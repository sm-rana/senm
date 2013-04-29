#include <stdio.h>
#include <Windows.h>
#include "toolkit.h"
#include "types.h"
#include "funcs.h"
#include "rngs.h" 
#include "rvgs.h" //random number generator
#include "DataSource.h"
#include "ScadaGenerator.h"


extern "C" { /* to link with epanet*/
	extern int d_update; //used to stop the EPS toolkit update demands 
	extern double* D;  //demands
}



ScadaGenerator::ScadaGenerator() : Provider() {

	n_chan = 0;
	chanlist = NULL;
	
}

ScadaGenerator::Err ScadaGenerator::init(
	const char* inpfile_path, 
	const char* dsn) {

	fprintf(stdout, "SCADA Data generator...\n\
					using inpfile: %s\n\
					using db connection: %s\n",
					inpfile_path, dsn);
	
	char tmp[MAX_FILE_PATH];
	strncpy_s(tmp, inpfile_path, MAX_FILE_PATH);
	/* Alloc mem for network components */
	if (ENopen(tmp, "", "")) return NET_FILE_ERR;

	/* alloc memory for demands*/
	ENgettimeparam(EN_DURATION, &dur);
	ENgettimeparam(EN_HYDSTEP, &hstep); // length of a hydraulic step
	ENgetcount(EN_NODECOUNT, &nnode);
	ENgetcount(EN_TANKCOUNT, &ntank);
	njunc = nnode - ntank;
	//nctrl = sizeof(ctrler)/sizeof(char*);  // number of controllable components

	fprintf(stdout, "The network has %d junctions, %d tanks/reserviors, \
					Simulation time span %d days, hydraulic time step is %d seconds.\n",
					njunc, ntank, dur/3600/24, hstep);

	if (hstep == 0 || njunc == 0) {
		fprintf(stderr, "Network error.\n");
		return NET_PROBLEM;
	}

	int nstep = dur/hstep;  // number of hydraulic steps
	demands = new double*[nstep + 1];
	int i = 0;
	for (i = 0; i < nstep + 1; ++i) {
		demands[i] = new double[njunc + 1];
	}

	ENopenH(); // Allocate mem for hyd simulation


	fprintf(stdout, "Run 1: Simulate hydraulics and obtain static demands\n");
	int istep = 0; // current hydraulic time step
	int problem;  // error code
	long step; /* epanet time until next simulation step*/
	long stime; /* epanet simulation clock in seconds */

	/* run simulation, get demand info */
	for (ENinitH(0), step=1; step>0; ENnextH(&step)) {
		// step - next hydraulic step length

		// fprintf(stdout, "Simulate hydraulics at Time %d second.\n",stime);
		problem = ENrunH(&stime);
		if (problem >= 100) {//severe errors in the original network
			fprintf(stdout, "Problem (%d): %s \n", problem, geterrmsg(problem));
			fprintf(stdout, "Quit!\n");
			return ORIG_NET_WONT_RUN;
		}

		/* pull demand info*/
		if (stime % hstep == 0) { //only get demand data in hstep-interval, (not on control step)
			for (i=1; i<=njunc; ++i) {
				demands[istep][i] = D[i];
			}
			istep++;
		}

	}
	
	fprintf(stdout, "Connecting to DB...\n");
	if ((ScadaGenerator::Err)connect(dsn)) return CANT_CONNECT_DB;

	/* Reset database schema */
	fprintf(stdout, "Building channels...\n");
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

	return OK;

}

ScadaGenerator::Err ScadaGenerator::make(
	SQL_TIMESTAMP_STRUCT* t1, int timespan, double cof) {

	SQLRETURN rc;	  /* db errors */
	int problem = 0;  /*epanet problems/errors */
		
	SelectStream(0);                  /* rngs: select the default stream */
	PutSeed(1); 

	int istep, i;  //timestep iterator
	int iround;  // round iterator
	double d, newd;  //temporary demand
	long stime = 0;  // epanet time
	int scada_time = 0;   // scada time shift
	long step; //time to next hyraulic simulation

	rc = SQLPrepareA(hStmt, (SQLCHAR*)"INSERT INTO Msmts (time, value, cid) \
					VALUES(DATE_ADD(?, INTERVAL ? SECOND), ?, ?);", SQL_NTS);
	if (rc) {
		handleDiagnosticRecord(hStmt, SQL_HANDLE_STMT, rc);
		if (rc == SQL_ERROR) return CANT_PREPARE_SQL;
	}

	SQLLEN cbdatetime2 = sizeof(SQL_TIMESTAMP_STRUCT);
	SQLLEN cbint = sizeof(int);
	SQLLEN cbfloat = sizeof(float);
	rc = SQLBindParameter(hStmt, 1, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, t1, 0, &cbdatetime2);
   

	/* simulating */

	for (iround = 0; iround < timespan/dur + 1; iround++ ) {
		// repeatedly run simulation until timespan is exhausted

		for (ENinitH(0), step=1, //  reset hydraulic engine
			istep=0, stime=0,				// reset hydraulic time/timestep iterator
			d_update=0;						// in epanet toolkit, disable auto demand update
			step>0; ENnextH(&step)) {

			scada_time = iround*dur + stime;
			rc = SQLBindParameter(hStmt, 2, SQL_PARAM_INPUT, 
				SQL_C_ULONG, SQL_INTEGER, 10, 0, &scada_time, 0, &cbint);

			fprintf(stdout, "Simulate hydraulics at Time %5.2f Hour. \n ",(float)stime/3600);

			if (stime % hstep == 0) { //only get demand data in hstep-interval, (not on control step)

				/* set water demands (Senario 1)*/
				for (i=1; i<=njunc; ++i) {
					d = demands[istep][i];
					if (d>0) { //water user
						newd = Lognormal2(d, cof*d);
						D[i] = newd; 
					}
				}
				istep++;
			}

			problem = ENrunH(&stime);
			if (problem >= 100) { // computational errors
				fprintf(stderr, "Problem (%d): %s \n", problem, geterrmsg(problem));
				return HYD_PROBLEM;
			}

			/* pull scada data*/
			if (scada_time <= timespan) //write db only when within timespan
			for (Channel* it=chanlist; it!=NULL; it=it->next) {

			/* go through all channels, write db*/
				rc = SQLBindParameter(hStmt, 4, SQL_PARAM_INPUT, 
				SQL_C_ULONG, SQL_INTEGER, 10, 0, &(it->key), 0, &cbint);

				float signal = 0;
				switch (it->mtype) {
					case Network::LSTATUS:
						ENgetlinkvalue(it->mindex, EN_STATUS, &signal);
						break;
					case Network::FLOW:
						ENgetlinkvalue(it->mindex, EN_FLOW, &signal);
						break;
					case Network::PRESSURE:
						ENgetnodevalue(it->mindex, EN_PRESSURE, &signal);
						break;
					case Network::HEAD:
						ENgetnodevalue(it->mindex, EN_HEAD, &signal);
						break;
				}

				rc = SQLBindParameter(hStmt, 3, SQL_PARAM_INPUT, 
				SQL_C_FLOAT, SQL_FLOAT, 15, 0, &signal, 0, &cbfloat);

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

}
