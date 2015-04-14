#include <stdio.h>
#include <tchar.h>
#include <strsafe.h>
#include <float.h>
#include <Windows.h>


#include "rvgs.h" //random number generator

#include "DataSource.h"
#include "ScadaGenerator.h"

#ifdef _DEBUG
#define debug 1
#else
#define debug 0
#endif

SG_ERR _loadNetAndDb(TCHAR* inpfile, char const* ini_path, ScadaGenerator** sg_out); 

SG_ERR SG_new(TCHAR* inp_path, char const* ini_path, int seed, double cov, ScadaGenerator** sg_out){
	SG_ERR err;

	err = _loadNetAndDb(inp_path, ini_path, sg_out);

	if (err != SG_OK) return err;

	(*sg_out)->rg = new Rvgs(seed);  // create random number generator
	(*sg_out)->mode = SG_IID_LN;
	(*sg_out)->_cov = cov;
	return SG_OK;
}

SG_ERR SG_new(TCHAR* inp_path, char const* ini_path, int seed, Varima* md, ScadaGenerator** sg_out){
	SG_ERR err;

	err = _loadNetAndDb(inp_path, ini_path, sg_out);

	if (err != SG_OK) return err;

	if (md == NULL || md->getN()!=(*sg_out)->nuser) {
		SG_ewi(err = SG_INVALID_VARIMA);
		return err;
	}

	int nstep = (*sg_out)->dur / (*sg_out)->hstep;
	(*sg_out)->rg = new Rvgs(seed);  // create random number generator
	if (nstep <= md->getLB()) {
		SG_ewi(err = SG_INIT_DATA_NOT_ENOUGH);
		return err;
	}
	md->initStateLogG( (*sg_out)->rg, (*sg_out)->_user, nstep);
	(*sg_out)->mode = SG_VARIMA;
	(*sg_out)->_md = md;

	return SG_OK;
}



SG_ERR _loadNetAndDb(TCHAR* inpfile, char const* ini_path, ScadaGenerator** sg_out) {
	USES_CONVERSION;
	SG_ERR err;

	ScadaGenerator* prec = (ScadaGenerator*)calloc(1, sizeof(ScadaGenerator));
	if (prec == NULL) return SG_MEM_NOT_ALLOCED;

	prec->rg = NULL;
	prec->_tpd = NULL;
	prec->_user = NULL;
	//prec->_net_array = NULL;
	//prec->_datastep = NULL;
	//	prec->_needsVis = 0;
	prec->demands = NULL; 
	prec->vdemands = NULL;
	prec->vpressures = NULL;
	prec->extObj = NULL;
	prec->extUpdate = NULL;

	/*Load the epanet network file */

	if (debug) _tprintf(TEXT("Creating SCADA generator using inpfile: %s and")
		TEXT(" data source configuration file: %s\n"), inpfile, A2T(ini_path));

	/* epanet toolkit needs ANSI string */
	char tmp[SG_MAX_FILE_PATH];
	int nFileNameLen = WideCharToMultiByte( CP_ACP, // ANSI Code Page
		0, // No special handling of unmapped chars
		inpfile, // wide-character string to be converted
		-1, NULL, 0, // No output buffer since we are calculating length
		NULL, NULL ); // Unrepresented char replacement - Use Default 
	if (nFileNameLen > SG_MAX_FILE_PATH) { //input file path too long
		SG_ewi(err = SG_FILE_PATH_TOO_LONG); goto END;}

	WideCharToMultiByte( CP_ACP, // ANSI Code Page
		0, // No special handling of unmapped chars
		inpfile, // wide-character string to be converted
		-1, tmp, SG_MAX_FILE_PATH, NULL, NULL );

	if (GetLastError() == ERROR_NO_UNICODE_TRANSLATION) {
		SG_ewi(err = SG_FILE_PATH_NOT_ANSI); goto END;} // epanet won't take non-ansi file names

	//1. open epanet (inp file)
	/* Alloc mem for network components */
	if (ENopen(tmp, "", "")) {
		SG_ewi(err = SG_NET_FILE_ERR); goto END;}

	/* alloc memory for demands*/
	ENgettimeparam(EN_DURATION, &prec->dur);
	ENgettimeparam(EN_HYDSTEP, &prec->hstep); // length of a hydraulic step
	ENgetcount(EN_NODECOUNT, &prec->nnode);
	ENgetcount(EN_TANKCOUNT, &prec->ntank);
	//int nlink;
	ENgetcount(EN_LINKCOUNT, &prec->nlink);
	prec->njunc = prec->nnode - prec->ntank;
	//nctrl = sizeof(ctrler)/sizeof(char*);  // number of controllable components

	if (debug) _tprintf(TEXT("The network has %d junctions, %d tanks/reserviors. ")
		TEXT("Simulation time span %d days, hydraulic time step is %d seconds.\n"),
		prec->njunc, prec->ntank, prec->dur/3600/24, prec->hstep);

	if (prec->hstep == 0 || prec->njunc == 0) {
		SG_ewi(err = SG_NET_PROBLEM); goto END;
	}

	int nstep = prec->dur/prec->hstep;  // number of hydraulic steps
	prec->demands = (double**)calloc(prec->njunc + 1, sizeof(double*));
	prec->vdemands = (double*)calloc(prec->nnode, sizeof(double));
	prec->vpressures = (double*)calloc(prec->nnode, sizeof(double));
	if (prec->demands == NULL || prec->vdemands == NULL || prec->vpressures == NULL) {
		SG_ewi(err = SG_MEM_NOT_ALLOCED); goto END;
	}

	int iJunc;
	for (iJunc = 0; iJunc <= prec->njunc; ++iJunc) {
		prec->demands[iJunc] = (double*)calloc(nstep + 1, sizeof(double));
		if (prec->demands[iJunc] == NULL) {
			SG_ewi(err = SG_MEM_NOT_ALLOCED); goto END;}
	}

	if (ENopenH()) { // Allocate mem for hyd simulation
		SG_ewi(err = SG_NET_INIT_ERR); goto END;}


	if (debug) _tprintf(TEXT("Run 1: Simulate hydraulics and obtain static demands\n"));
	int istep = 0; // current hydraulic time step
	int problem;  // error code
	long step; /* epanet time until next simulation step*/
	long stime; /* epanet simulation clock in seconds */

	prec->ntstep = 0;
	prec->maxd = -DBL_MAX;
	prec->mind = DBL_MAX;

	/* test-run simulation, get demand info */
	for (ENinitH(0), step=1; step>0; ENnextH(&step)) {
		// step - next hydraulic step length

		// fprintf(stdout, "Simulate hydraulics at Time %d second.\n",stime);
		problem = ENrunH(&stime);
		if (problem >= 100) {//severe errors in the original network
			if (debug) _tprintf(TEXT("Problem (%d): %s \n"), problem, geterrmsg(problem));
			if (debug) _tprintf(TEXT("Quit!\n"));
			SG_ewi(SG_ORIG_NET_WONT_RUN);
			goto END;
		}

		/* pull demand info*/
		if (stime >= istep*prec->hstep) { 
			//only get demand data in hstep-interval, (not on control step)
			for (iJunc=1; iJunc<=prec->njunc; ++iJunc) {
				prec->demands[iJunc][istep] = D[iJunc];  
				if (D[iJunc] > prec->maxd ) prec->maxd = D[iJunc];
				if (D[iJunc] < prec->mind ) prec->mind = D[iJunc];
			}
			istep++;
		}
		prec->ntstep++;
	}
    ENcloseH();

	// 2. prepare db connection
	if (ini_path) {

		if (debug) _tprintf(TEXT("Connecting to DB and building channels...\n"));
		if (DataSource::New(ini_path, &prec->ds) != DataSource::OK ) {
			SG_ewi(err = SG_CANT_EST_DB); goto END;}
		if (prec->ds->n_prov > 1) {
			SG_ewi(err = SG_MORE_THAN_ONE_PROVIDER); goto END;}
		if (prec->ds->lsProv->resetFactTab() != DataSource::OK) {//clean existing data
			SG_ewi(err = SG_CANT_EST_DB); goto END;}

		/* Verify channel net_id and find out index */
		int index;  // network index
		for (Channel* it = prec->ds->lsChan; it != NULL; it = it->next) {
			if (ENgetlinkindex(it->name, &index)) { // can't find the link
				if (ENgetnodeindex(it->name, &index)) { // can't find the node
					SG_ewi(err = SG_UNKNOWN_CHANNEL); goto END;}
			}
			it->mindex = index;
		}
	} else {
		if (debug) _tprintf(TEXT("Dry-run: no database connection"));
		prec->ds = NULL;
	}

	// 3. prepare water demand storage, demand model interface
	prec->nuser = 0;
	/* set nuser*/
	for (iJunc=1; iJunc<=prec->njunc; ++iJunc) {
		double tp = 0; // total demand for a junction in the duration
		for (istep=0; istep<nstep; ++istep) 
			tp += prec->demands[iJunc][istep];
		if (tp!=0) {  // is a water user
			++prec->nuser;
		} 
	}

	prec->_user = (double**)calloc(prec->njunc+1, sizeof(double*)); 
	int i; // water demands iter

	//alloc njunc+1, but only nuser+1 will be used at most
	int j;  // water user iterator
	for (i=1, j=0 ; i<=prec->njunc; ++i) {
		double tp = 0; // total demand for a junction in the duration
		for (istep=0; istep<nstep; ++istep) tp += prec->demands[i][istep];
		if (tp!=0) {  // is a water user
			prec->_user[j++] = prec->demands[i];
		} else {
			free(prec->demands[i]); // release mem
			prec->demands[i] = NULL;
		}
	}

	prec->_tpd = (double*)calloc(prec->nuser, sizeof(double)); //

	*sg_out = prec;
	return SG_OK;

END: //clean up when an error is raised
	if (prec->demands) {
		for (iJunc = 1; iJunc<=prec->njunc; ++iJunc)
			free(prec->demands[iJunc]);
	}
	free(prec->demands);
	free(prec->vdemands);
	free(prec->vpressures);
	free(prec->_user);
	free(prec->_tpd);
	if (prec->rg) delete prec->rg;
	free(prec);
	return err;
}


void SG_ewi(SG_ERR err) {
	TCHAR errorTxts[SG_DUMMY_LAST][MAX_ERR_STRING_SIZE-32] = {
		TEXT("OK"),

		TEXT("Error in reading EPANET INP file."),
		TEXT("No EPS simulation or no juntions in the EPANET network."),
		TEXT("Could not initialize EPANET network (mem allocation error?)."),
		TEXT("Could not run EPANET hydraulic simulation."),
		TEXT("Could not connect to the database or database broken."),
		TEXT("Data source contains more than one provider."),

		TEXT("Could not initialize scada database."),
		TEXT("Channel information does not match EPANET network"),
		TEXT("Errors in hydraulic simulation."),
		TEXT("Could not write to the scada database."),

		TEXT("Given INP file name is too long."),
		TEXT("Given INP file path/name is not in ASCII"),
		TEXT(""),
		TEXT(""),

		TEXT(""),
		TEXT("Memory not allocated."),
		TEXT("Data used for initliazation is not long enough for the demand model."),
	};

	TCHAR errTxt[MAX_ERR_STRING_SIZE];
	_stprintf(errTxt, TEXT("<GEN>ScadaGenerator: %s"), errorTxts[err]);

	ewi(errTxt);
}






SG_ERR SG_make(ScadaGenerator *sg, Tstamp* t0_in, int timespan) {

	int problem = 0;  /*epanet problems/errors */
	SG_ERR err;


	int istep, i;  //timestep iterator
	int iround;  // round iterator
	double d, newd;  //temporary demand
	long stime = 0;  // epanet time
	int scada_time_shift = 0;   // scada time shift
	long step; //time to next hyraulic simulation
	int itstep;
	sg->elapTime = 0;
	sg->t0 = CTime(t0_in->year, t0_in->month, t0_in->day, t0_in->hour, 
		t0_in->minute, t0_in->second, 0);

	if (sg->ds) {// has db connection
		if (sg->ds->lsProv->initInsertScada(t0_in) != Provider::OK) {
			SG_ewi(err = SG_SCADA_INIT_ERR); return err;
		}
	}

	// repeatedly run simulation until timespan is met
	for (iround = 0; iround < timespan/sg->dur + 1; iround++ ) 
		for (ENopenH(), ENinitH(0), step=1, 
			/*  reset the hydraulic engine, here we assume the original network 
			has been calibrated so that the tank levels at the end of simulation
			match those at the beginning of the simulation */
			istep=0, stime=0,	// reset hydraulic time/timestep iterator
			itstep = 0,         // reset total step counter
			d_update=0;			// in epanet toolkit, disable auto demand update
	step>0; ENnextH(&step), ++itstep) {

		long tempt;
		ENgettimeparam(EN_HTIME, &tempt);
		sg->elapTime = iround*sg->dur + tempt;

		if (stime + step == sg->dur) break; //don't compute last snapshot-will overlap with next round

		if (stime + step >= (istep)*sg->hstep) { 
			//only set demand data in hstep-interval, (not at control step)

			switch (sg->mode) {
			case SG_IID_LN:
				/* set water demands (Senario 1)*/
				for (i=1; i<=sg->njunc; ++i) {
					d = sg->demands[i][istep];
					if (d>0) { //water user
						newd = sg->rg->Lognormal2(d, sg->_cov*d);
						D[i] = newd; 
					}			
				}
				break;
			case SG_VARIMA:
				sg->_md->generateExp(sg->_tpd);
				int iuser;
				for (i=1, iuser = 0; i<=sg->njunc; ++i)
					if (sg->demands[i])  D[i] = sg->_tpd[iuser++];
				break;
			}
			istep++;
		} 

		// run hydraulic simulation
		problem = ENrunH(&stime);
		scada_time_shift = iround*sg->dur + stime;

		if (debug) _ftprintf(stderr, TEXT(
			"Hyd-sim @ %5.2f Hr (+ %d * %d s) \n "),
			(float)stime/3600, iround, sg->dur);

		if (problem >= 100) { // computational errors
			if (debug) _ftprintf(stderr, 
				TEXT("Problem (%d): %s \n"), problem, geterrmsg(problem));
			SG_ewi(err = SG_HYD_PROBLEM);
			return err;
		}


		/* pull scada data*/
		if (scada_time_shift <= timespan && sg->ds) {//write db only when within timespan
			Channel* it; 
			for (it=sg->ds->lsChan; it!=NULL; it=it->next) {

				float signal = 0;
				switch (it->type) {
				case Channel::C:
					ENgetlinkvalue(it->mindex, EN_STATUS, &signal);
					break;
				case Channel::Q: case Channel::F: //pipe and pump flow rate
					ENgetlinkvalue(it->mindex, EN_FLOW, &signal);
					break;
				case Channel::P: case Channel::B: case Channel::A:
					ENgetnodevalue(it->mindex, EN_PRESSURE, &signal);
					break;
				case Channel::L:
					ENgetnodevalue(it->mindex, EN_HEAD, &signal);
					break;
				case Channel::D:
					ENgetnodevalue(it->mindex, EN_DEMAND, &signal);
					break;
				case Channel::V:
					ENgetlinkvalue(it->mindex, EN_SETTING, &signal);
					break;
				default:
					continue;
				}
				if (sg->ds->lsProv->insertScada(&scada_time_shift, &it->key, &signal) != Provider::OK) {
					SG_ewi(err = SG_CANT_FILL_SCADA); return err; }
			}
		}

		if (sg->extObj && sg->extUpdate) {//extra function call, can be used for visualization
			for (int ii=1; ii<=sg->nnode; ++ii) {
				float tpp1, tpp2;
				ENgetnodevalue(ii, EN_HEAD, &tpp1);
				ENgetnodevalue(ii, EN_ELEVATION, &tpp2);
				sg->vpressures[ii-1] = tpp1-tpp2;
				sg->vdemands[ii-1] = ii>sg->njunc?0:D[ii];
			}

			(*sg->extUpdate)(sg->extObj, sg);
		}


	}

	ENcloseH(); // free mem used for hyd simulation
	return SG_OK;
}



SG_ERR SG_delete(ScadaGenerator *sg) {
	ENclose();  //free mem used for net components

	free(sg->_user);
	free(sg->_tpd);
	free(sg->vdemands);
	free(sg->vpressures);
	for (int i=0;i<sg->njunc;++i)
		free(sg->demands[i]);
	free(sg->demands);

	delete sg->ds;

	delete sg->rg;

	return SG_OK;
}
