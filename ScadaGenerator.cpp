#include <stdio.h>
#include <Windows.h>
#include "toolkit.h"
#include "types.h"
#include "funcs.h"
#include "rvgs.h" //random number generator

#include "DataSource.h"
#include "ScadaGenerator.h"
#include <tchar.h>
#include <strsafe.h>

extern "C" { /* to link with epanet*/
	//prevent the EPS func from updating demands, see epanet/hydraul.c
	extern int d_update; 

	//gain access to water demands directly
	extern double* D;  
	extern Snode* Node;
	extern Slink* Link;
}

ScadaGenerator::ScadaGenerator() : Provider() {

	ntstep = 0;
	n_chan = 0;
	chanlist = NULL;
	rg = NULL;
	nuser = 0;
	_tpd = NULL;
	_user = NULL;
	_net_array = NULL;
	_datastep = NULL;
	_needsVis = 0;
	minx=miny=mind=minz=BIG;
	maxx=maxy=maxd=maxz=-BIG;
	write_db = 1;
	
}
ScadaGenerator::Err ScadaGenerator::init(
	LPCTSTR inpfile_path, 
	LPCTSTR dsn, int seed) {
	return init(inpfile_path, dsn, seed, NULL);
}

ScadaGenerator::Err ScadaGenerator::init(
	LPCTSTR inpfile_path, 
	LPCTSTR dsn, int seed, vtkPolyData* nettop[2]) {

	/* initialize the database and load the epanet network file */

	_tprintf(TEXT("SCADA Data generator...\n\
					using inpfile: %s\n\
					using db connection: %s\n"),
					inpfile_path, dsn);
	
	/* epanet toolkit needs ANSI string */
	char tmp[MAX_FILE_PATH];
	int nFileNameLen = WideCharToMultiByte( CP_ACP, // ANSI Code Page
		0, // No special handling of unmapped chars
		inpfile_path, // wide-character string to be converted
		-1,
		NULL, 0, // No output buffer since we are calculating length
		NULL, NULL ); // Unrepresented char replacement - Use Default 
	if (nFileNameLen > MAX_FILE_PATH)  //input file path too long
		return FILE_PATH_TOO_LONG;
	WideCharToMultiByte( CP_ACP, // ANSI Code Page
		0, // No special handling of unmapped chars
		inpfile_path, // wide-character string to be converted
		-1,
		tmp, 
		MAX_FILE_PATH,
		NULL, NULL );

	if (GetLastError() == ERROR_NO_UNICODE_TRANSLATION)
		return FILE_PATH_NOT_ANSI;  // epanet won't take non-ansi file names

	/* Alloc mem for network components */
	if (ENopen(tmp, "", "")) return NET_FILE_ERR;

	/* alloc memory for demands*/
	ENgettimeparam(EN_DURATION, &dur);
	ENgettimeparam(EN_HYDSTEP, &hstep); // length of a hydraulic step
	ENgetcount(EN_NODECOUNT, &nnode);
	ENgetcount(EN_TANKCOUNT, &ntank);
	int nlink;
	ENgetcount(EN_LINKCOUNT, &nlink);
	njunc = nnode - ntank;
	//nctrl = sizeof(ctrler)/sizeof(char*);  // number of controllable components

	// set coordinates
	int i; 
	if (nettop!=NULL) {
				
		double coords[3];

		vtkPoints* pts3 = vtkPoints::New();
		vtkPoints* pts2 = vtkPoints::New();
		vtkCellArray* ca = vtkCellArray::New();  // cells for nodes
		vtkCellArray* cl = vtkCellArray::New();  // cells for links
		
		_netpoly[0] = vtkPolyData::New();
		_netpoly[1] = vtkPolyData::New(); //flat
		vtkPolyData* tpn1 = vtkPolyData::New();

		for (i=1; i<=nnode; ++i) {
			coords[0] = Node[i].x;
			if (coords[0] > maxx ) maxx = coords[0]; //update bands
			if (coords[0] < minx ) minx = coords[0];
			coords[1] = Node[i].y;
			if (coords[1] > maxy ) maxy = coords[1];
			if (coords[1] < miny ) miny = coords[1];
			coords[2] = Node[i].El;
			if (coords[2] > maxz ) maxz = coords[2];
			if (coords[2] < minz ) minz = coords[2];

			pts3->InsertPoint(i-1, coords);
			vtkIdType icell = i-1;
			ca->InsertNextCell(1, &icell);

			// flat map
			coords[2] = 0;
			pts2->InsertPoint(i-1, coords);

		}

		_netpoly[0]->SetPoints(pts3);
		_netpoly[1]->SetPoints(pts2); 
		_netpoly[0]->SetVerts(ca);
		_netpoly[1]->SetVerts(ca);

		for (i=1; i<nlink; ++i) {
			vtkIdType jcell[2];
			jcell[0] = Link[i].N1 - 1; // to vis index
			jcell[1] = Link[i].N2 - 1;
			cl->InsertNextCell(2, jcell);
		}

		_netpoly[0]->SetLines(cl);
		_netpoly[1]->SetLines(cl);


		vtkTransform* bz = vtkTransform::New();
			bz->Identity();
			bz->Scale(1, 1, 1 /*(maxx-minx)/3/(maxz-minz)*/);
		vtkTransformPolyDataFilter* tpf = vtkTransformPolyDataFilter::New();
			tpf->SetTransform(bz);
			tpf->AddInput(_netpoly[0]);
			tpf->Update();

	
		
		nettop[0] = tpf->GetOutput();
		nettop[1] = _netpoly[1];

	}
	
	_tprintf(TEXT("The network has %d junctions, %d tanks/reserviors, \
			Simulation time span %d days, hydraulic time step is %d seconds.\n"),
			njunc, ntank, dur/3600/24, hstep);

	if (hstep == 0 || njunc == 0) {
		_tprintf(TEXT("Network error.\n"));
		return NET_PROBLEM;
	}

	int nstep = dur/hstep;  // number of hydraulic steps
	demands = new double*[njunc + 1];

	for (i = 0; i < njunc + 1; ++i) {
		demands[i] = new double[nstep + 1];
	}

	ENopenH(); // Allocate mem for hyd simulation


	_tprintf(TEXT("Run 1: Simulate hydraulics and obtain static demands\n"));
	int istep = 0; // current hydraulic time step
	int problem;  // error code
	long step; /* epanet time until next simulation step*/
	long stime; /* epanet simulation clock in seconds */

	ntstep = 0;
	/* run simulation, get demand info */
	for (ENinitH(0), step=1; step>0; ENnextH(&step)) {
		// step - next hydraulic step length

		// fprintf(stdout, "Simulate hydraulics at Time %d second.\n",stime);
		problem = ENrunH(&stime);
		if (problem >= 100) {//severe errors in the original network
			_tprintf(TEXT("Problem (%d): %s \n"), problem, geterrmsg(problem));
			_tprintf(TEXT("Quit!\n"));
			return ORIG_NET_WONT_RUN;
		}

		/* pull demand info*/
		if (stime >= istep*hstep) { //only get demand data in hstep-interval, (not on control step)
			for (i=1; i<=njunc; ++i) {
					demands[i][istep] = D[i];  
					if (D[i] > maxd ) maxd = D[i];
					if (D[i] < mind ) mind = D[i];
			}
			istep++;
			
		}

		ntstep++;

	}

	if (dsn) {

		_tprintf(TEXT("Connecting to DB...\n"));
		if ((ScadaGenerator::Err)connect(dsn)) return CANT_CONNECT_DB;

		/* Reset database schema */
		_ftprintf(stdout, TEXT("Building channels...\n"));
		if ((ScadaGenerator::Err)loadChannels(&chanlist, &n_chan)) 
			return CHAN_TABLE_ERR;

		/* Verify channel net_id and find out index */
		int index;  // network index
		for (Channel* it = chanlist; it != NULL; it = it->next) {
			if (ENgetlinkindex(it->name, &index)) { // can't find the link
				if (ENgetnodeindex(it->name, &index)) { // can't find the node
					return UNKNOWN_CHANNEL;
				}
			}
			it->mindex = index;
		}

	} else {
		write_db = 0;
	}

	/* set nuser*/
	for (i=1; i<=njunc; ++i) {
		double tp = 0; // total demand for a junction in the duration
		for (istep=0; istep<nstep; ++istep) tp += demands[i][istep];
		if (tp!=0) {  // is a water user
			++nuser;
		} 
	}

	rg = new Rvgs(seed);  // create random number generator

	return OK;


}

ScadaGenerator::Err ScadaGenerator::make(
	Tstamp* t1, int timespan, Varima* md) {
	
	return make(t1, timespan, md, NULL, NULL, NULL);
}


ScadaGenerator::Err ScadaGenerator::make(
	Tstamp* t1, int timespan, 
	Varima* md, vtkPolyData** arr, volatile LONG* datastep, int* dispstep) {

	mode = TEMPORAL_SPATIAL_CORRELATED;

	if (md == NULL ) return INVALID_VARIMA;
	if (arr != NULL && datastep!= NULL && dispstep != NULL) {
		_needsVis = 1;
		_net_array = arr;
		_datastep = datastep;
		_dispstep = dispstep;
		vtkCylinderSource* cyl0 = vtkCylinderSource::New();
		cyl0->SetHeight(5);
		cyl0->SetRadius(1);
		cyl0->SetCenter(0, -2.5, 0);
		vtkTransform* ytoz = vtkTransform::New();
			ytoz->Identity();
			ytoz->RotateX(90);
		cyl = vtkTransformFilter::New();
			cyl->SetTransform(ytoz);
			cyl->AddInputConnection(cyl0->GetOutputPort());
			cyl->Update();

		zlineP = vtkPolyData::New();
		zlineD = vtkPolyData::New();
		vtkPoints* ptsP = vtkPoints::New();
		vtkPoints* ptsD = vtkPoints::New();
		ptsP->InsertPoint(0, 0, 0, 0);
		ptsD->InsertPoint(0, 0, 0, 0);
		ptsP->InsertPoint(1, 0, 0, (maxx-minx+maxy-miny)/4/(maxz-minz)); //set scale factor
		ptsD->InsertPoint(1, 0, 0, (maxx-minx+maxy-miny)/4/(maxd-mind));
		vtkCellArray* cl = vtkCellArray::New();
		vtkIdType jcell[] = {0, 1};
		cl->InsertNextCell(2, jcell);
		zlineP->SetPoints(ptsP);
		zlineD->SetPoints(ptsD);
		zlineP->SetLines(cl);
		zlineD->SetLines(cl);

	}
	
	int i, istep; // water demands iter
	int nstep = dur/hstep;

	_user = (double**)calloc(njunc+1, sizeof(double*)); 
	//alloc njunc+1, but only nuser+1 will be used at max

	int j;  // water user iterator
	for (i=1, j=0 ; i<=njunc; ++i) {
		double tp = 0; // total demand for a junction in the duration
		for (istep=0; istep<nstep; ++istep) tp += demands[i][istep];
		if (tp!=0) {  // is a water user
			_user[j++] = demands[i];
		} else {
			delete[] demands[i]; // release mem
			demands[i] = NULL;
		}
	}

	_tpd = (double*)calloc(nuser, sizeof(double)); //

	if (md->getN() != nuser) return VARIMA_N_NOT_MATCH;
	// fill _use matrix, so junctions w/o demands are not included in varima	
	
	_md = md;
	_md->initStateLogG(rg, _user, nstep); // log(demands) to init the arima model
	
	return make(t1, timespan);	

}
ScadaGenerator::Err ScadaGenerator::make(
	Tstamp* t1, int timespan, double cov) {

		mode = INDEPENDENT;
		_cov = cov;

		return make(t1, timespan);

}
ScadaGenerator::Err ScadaGenerator::make(
	Tstamp* t1, int timespan) {

	SQLRETURN rc;	  /* db errors */
	int problem = 0;  /*epanet problems/errors */
	

	int istep, i;  //timestep iterator
	int iround;  // round iterator
	double d, newd;  //temporary demand
	long stime = 0;  // epanet time
	int scada_time = 0;   // scada time shift
	long step; //time to next hyraulic simulation
	
	SQLTCHAR sqltext[MAX_SQL_LEN];
	_stprintf(sqltext, TEXT("INSERT INTO %s (time, value, cid) \
					VALUES(DATE_ADD(?, INTERVAL ? SECOND), ?, ?);"), dat_tab);

	rc = SQLPrepare(hStmt, sqltext, SQL_NTS);
	if (rc) {
		handleDiagnosticRecord(hStmt, SQL_HANDLE_STMT, rc);
		if (rc == SQL_ERROR) return CANT_PREPARE_SQL;
	}

	rc = SQLBindParameter(hStmt, 1, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, t1, 0, NULL);
   
	/* simulating */
	
	int itstep;
	for (iround = 0; iround < timespan/dur + 1; iround++ ) {
		// repeatedly run simulation until timespan is met

		for (ENinitH(0), step=1, 
			/*  reset the hydraulic engine, here we assume the original network 
			has been calibrated so that the tank levels at the end of simulation
			match those at the beginning of the simulation */
								 
			istep=0, stime=0,	// reset hydraulic time/timestep iterator
			itstep = 0,         // reset total step counter
			d_update=0;			// in epanet toolkit, disable auto demand update
			step>0; ENnextH(&step), ++itstep) {

			if (stime == dur) break; //don't compute last snapshot
			scada_time = iround*dur + stime;
			rc = SQLBindParameter(hStmt, 2, SQL_PARAM_INPUT, 
				SQL_C_ULONG, SQL_INTEGER, 10, 0, &scada_time, 0, NULL);

			_ftprintf(stdout, TEXT("Run hydraulics at Time %5.2f Hr. \n "),(float)stime/3600);

			if (stime + step >= (istep)*hstep) { 
				//only set demand data in hstep-interval, (not at control step)

				switch (mode) {
				case INDEPENDENT:
					/* set water demands (Senario 1)*/
					for (i=1; i<=njunc; ++i) {
						d = demands[i][istep];
						if (d>0) { //water user

							newd = rg->Lognormal2(d, _cov*d);
							D[i] = newd; 
						}			
					}
					break;
				case TEMPORAL_SPATIAL_CORRELATED:
					_md->generateExp(_tpd);
					int iuser;
					if (d_update == 0) //epanet demand update is disabled
					for (i=1, iuser = 0; i<=njunc; ++i)
						if (demands[i]) {
							D[i] = _tpd[iuser++];
						//	if (_needsVis) dmds->InsertTuple1(i-1, D[i]);
						} else {
						//	if (_needsVis) dmds->InsertTuple1(i-1, 0);
						}
					
					break;
					
				}
				istep++;
			} 

			problem = ENrunH(&stime);
			if (problem >= 100) { // computational errors
				_ftprintf(stderr, TEXT("Problem (%d): %s \n"), problem, geterrmsg(problem));
				return HYD_PROBLEM;
			}

			if (_needsVis) { 
				//vis, to run this part, _net_array must be allocated
				vtkDoubleArray* pressures = vtkDoubleArray::New();
				vtkDoubleArray* dmds = vtkDoubleArray::New();
				for (int ii=1; ii<=nnode; ++ii) {
					float tpp1, tpp2;
					ENgetnodevalue(ii, EN_HEAD, &tpp1);
					ENgetnodevalue(ii, EN_ELEVATION, &tpp2);
					
					pressures->InsertTuple1(ii-1, tpp1-tpp2);
					
					dmds->InsertTuple1(ii-1, ii>njunc? 0 :D[ii]);
					 
				}

				vtkPolyData* tmpoly1 = vtkPolyData::New(); 
				vtkPolyData* tmpoly2 = vtkPolyData::New();
				tmpoly1->DeepCopy(_netpoly[0]);
				tmpoly1->GetPointData()->SetScalars(pressures);
				pressures->Delete();
				tmpoly2->DeepCopy(_netpoly[1]);
				tmpoly2->GetPointData()->SetScalars(dmds);
				dmds->Delete();

				vtkGlyph3D* gly = vtkGlyph3D::New();
				vtkGlyph3D* gly2 = vtkGlyph3D::New();
				gly->SetInput(tmpoly1);
				gly2->SetInput(tmpoly2);
				//gly->SetSourceConnection(cyl->GetOutputPort());
				//gly2->SetSourceConnection(cyl->GetOutputPort());
				gly->SetSource(zlineP);
				gly2->SetSource(zlineD);
				gly->SetScaleModeToScaleByScalar();
				gly2->SetScaleModeToScaleByScalar();
				//gly->Update();
				//gly2->Update();

				vtkTubeFilter* tubP = vtkTubeFilter::New();
				vtkTubeFilter* tubD = vtkTubeFilter::New();
				tubP->SetInputConnection(gly->GetOutputPort());
				tubP->SetNumberOfSides(10);
				tubP->SetRadius((maxx-minx)/400);
				tubD->SetInputConnection(gly2->GetOutputPort());
				tubD->SetNumberOfSides(10);
				tubD->SetRadius((maxx-minx)/100);
				tubP->SetCapping(1); 
				tubD->SetCapping(1); 
				
				tubP->Update();
				tubD->Update();

				_net_array[(*_datastep+1)*2] = tubP->GetOutput();
				_net_array[(*_datastep+1)*2+1] = tubD->GetOutput();

				// display time step
				long tempt;
				ENgettimeparam(EN_HTIME, &tempt);
				_dispstep[*_datastep+1] = iround*dur + tempt;

				InterlockedAdd(_datastep, 1);
			}

			/* pull scada data*/
			if (scada_time <= timespan && write_db) //write db only when within timespan
			for (Channel* it=chanlist; it!=NULL; it=it->next) {

			/* go through all channels, write db*/
				rc = SQLBindParameter(hStmt, 4, SQL_PARAM_INPUT, 
				SQL_C_ULONG, SQL_INTEGER, 10, 0, &(it->key), 0, NULL);

				float signal = 0;
				switch (it->type) {
				case Channel::C:
						ENgetlinkvalue(it->mindex, EN_STATUS, &signal);
						break;
				case Channel::Q:
						ENgetlinkvalue(it->mindex, EN_FLOW, &signal);
						break;
				case Channel::P:
						ENgetnodevalue(it->mindex, EN_PRESSURE, &signal);
						break;
				case Channel::L:
						ENgetnodevalue(it->mindex, EN_HEAD, &signal);
						break;
				}

				rc = SQLBindParameter(hStmt, 3, SQL_PARAM_INPUT, 
				SQL_C_FLOAT, SQL_FLOAT, 15, 0, &signal, 0, NULL);

				rc = SQLExecute(hStmt);
				if (rc) {
					handleDiagnosticRecord(hStmt, SQL_HANDLE_STMT, rc);
					if (rc == SQL_ERROR) return CANT_FILL_SCADA;
				}
		
			}

		}
	}

	return OK;
}


ScadaGenerator::~ScadaGenerator() {
	ENcloseH(); // free mem used for hyd simulation
	ENclose();  //free mem used for net components

	free(_user);
	free(_tpd);

	delete(rg);
	delete[] demands;

}
