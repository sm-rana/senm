#include <stdio.h>
#include "toolkit.h"
#include "types.h"
#include "funcs.h"
#include "rngs.h" 
#include "rvgs.h" //random number generator
#include "DataSource.h"

/* this routine fill the database with generated SCADA data 
from a network and pre-defined temporal-spatial correlation
of water demands */

extern "C" {
	extern int d_update; //used to stop the EPS toolkit update demands 
	extern double* D;  //demands
}

int main1 ()  {

	char* inp = "Data\\c-town_true_network.inp";
	
	char* ctrler[] = {"PU1", "PU2", "V2", "PU4", "PU5", "PU6", "PU7", 
		"PU8", "PU10", "PU11"};
	double cof = 0.05;  //coefficient of variance set for water demands

	CONST WCHAR* dsn = 
		L"DSN=jcx;\
		 DESCRIPTION={Jinduan's Senm Database};\
		 SERVER=localhost;\
		 UID=rtx_db_agent;\
		 PWD=rtx_db_agent;\
		 DATABASE=rtx_demo_db;\
		 PORT=3306;\
		 FOUND_ROWS=1";

	Provider* db = new Provider;

	fprintf(stdout, "SCADA Data generator...\n\
					using inpfile: %s\n\
					using db connection: %s\n",
					inp, dsn);
	ENopen(inp, "", "");

	long step(1), stime(0), dur(0), hstep(0);
	int nnode, ntank, njunc, nctrl;


	/* alloc memory for demands*/
	ENgettimeparam(EN_DURATION, &dur);
	ENgettimeparam(EN_HYDSTEP, &hstep); // length of a hydraulic step
	ENgetcount(EN_NODECOUNT, &nnode);
	ENgetcount(EN_TANKCOUNT, &ntank);
	njunc = nnode - ntank;
	nctrl = sizeof(ctrler)/sizeof(char*);  // number of controllable components

	fprintf(stdout, "The network has %d junctions, %d tanks/reserviors, \
					Simulation time span %d days, hydraulic time step is %d seconds.\n",
					njunc, ntank, dur/3600/24, hstep);

	if (hstep == 0 || njunc == 0) {
		fprintf(stderr, "Network error.\n");
		return 1;
	}

	int nstep = dur/hstep;  // number of hydraulic steps
	double** demands = new double*[nstep + 1];
	char**	controls = new char*[nstep + 1]; 
	int i = 0;
	for (i = 0; i < nstep + 1; ++i) {
		demands[i] = new double[njunc + 1];
		controls[i] = new char[nctrl];
	}


	fprintf(stdout, "Run 1: Simulate hydraulics and obtain static demands\n");
	int istep = 0; // current hydraulic time step
	int problem;  // error code

	ENopenH();
	/* run simulation, get demand info */
	for (ENinitH(0), step=1; step>0; ENnextH(&step)) {
		// step - next hydraulic step length

		fprintf(stdout, "Simulate hydraulics at Time %d second.\n",stime);


		problem = ENrunH(&stime);
		if (problem >= 100) {//severe errors in the original network
			fprintf(stdout, "Problem (%d): %s \n", problem, geterrmsg(problem));
			fprintf(stdout, "Quit!\n");
			goto cleanup;
		}

		/* pull demand info*/
		if (stime % hstep == 0) { //only get demand data in hstep-interval, (not on control step)
			for (i=1; i<=njunc; ++i) {
				demands[istep][i] = D[i];
			}
			istep++;
		}

	}
	// ENcloseH();

	fprintf(stdout, "Run 2: Simulate hydraulics with stochastic demands\n");
	/* run simulation with randomized demands */
	// Possible randomization scenarios:
	// 1. independent log-normal rv for each node at each step
	// 2. temporally correlated demands (e.g, log-AR(2))
	// 3. spatially correlated demands (
	// 4. temporally and spatially correlated

	SelectStream(0);                  /* rngs: select the default stream */
	PutSeed(1); 

	double d, newd;  //temporary demand
	for (ENinitH(0), step=1, //  reset hydraulic engine
		istep=0, stime=0,				// reset hydraulic time/timestep iterator
		d_update=0;						// in epanet toolkit, disable auto demand update
		step>0; ENnextH(&step)) {

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
				fprintf(stdout, "Problem (%d): %s \n", problem, geterrmsg(problem));
			}

			/* pull control info*/
			for (i=0; i<nctrl; ++i) {
				int index;
				float status;
				ENgetlinkindex(ctrler[i], &index);
				ENgetlinkvalue(index, EN_STATUS, &status);
				controls[istep][i] = (char)status;
			}

	}


cleanup:
	ENcloseH();
	ENclose();
	for (int i = 0; i < nstep + 1; ++i) {
		delete demands[i];
	} 
	delete demands;
	delete db;



}