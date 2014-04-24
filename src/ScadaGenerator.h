/* SenM - Sensor-driven hydraulic Model of water distribution systems
Copyright (C) 2013  Jinduan Chen <jinduan.uc@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License v2
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html)
for more details.   */

#pragma once

#include <tchar.h>
#include <process.h>
#include "DataSource.h"
#include "Varima.h"
#include "SenmCoreIncs.h"
#include "toolkit.h"
#include "types.h"
#include "funcs.h"

#ifdef _DEBUG
#define debug 1
#else
#define debug 0
#endif

#define SG_MAX_FILE_PATH  255 ///> max number of characters in a file path
#define SG_MAX_ERR_TEXT 512  ///> max size of error text

// Get access to some EPANET data structures
extern "C" { /* to link with epanet*/
	//prevent the EPS func from updating demands, see epanet/hydraul.c
	extern int d_update; 

	//gain access to water demands directly
	extern double* D;  
	extern Snode* Node;
	extern Slink* Link;
}


enum SG_Mode { 
	/// independent log-normal rv for each node at each step
	SG_IID_LN,
	/// vector time series (i.e., temporally and spatially correlated)
	SG_VARIMA, 
};
/// Generate synthetic SCADA data using EPANET
/*A ScadaGenerator populates the database with generated SCADA data 
from a network file and a demand model (can be set as a iid mutlti-dimensional
normal distribution model or a seasonal vectorized time series model) .
*/

struct ScadaGenerator{
	/** Types of water demand model*/

	//TCHAR* _inp; ///> input file name
	SG_Mode mode;  ///> generation mode 
	DataSource* ds;  ///> db data source, if NULL, dry-run without db

	long dur; /* epanet simulation duration in seconds */
	long hstep; /* epaent hydraulic time step length */
	int nnode, ntank, njunc, nlink; /* counts of nodes, tanks, junctions*/
	int nuser;  //no. of water users, water user must have demand >0 
	//at at least one time step
	int ntstep; //total time steps/snapshots, including control steps
	double mind, maxd; //minimum and maximum water demand 

    CTime t0; // simulated scada record starting time
	int elapTime; // elapsed time since beginning of the simulation


	double** demands;  //storage for mean demands at all hyd steps
	//demand[i][istep] is the ith nodal demand at timestep istep

	/// coefficient of variance. for SG_IID_LN mode only
	double _cov;
	Rvgs* rg;  // random num generator, 

	// varima model, for VARIMA mode only
	Varima* _md;

	//map from water user id to water demands (junction id) 
	//pointers to water users'
	double** _user;  

	// temporary water demand storage for synthetic water demands
	double* _tpd; 

	// snapshot of demands and presssures - to be used for visualization
	double* vdemands;
	double* vpressures;

	// function hook, will be executed during each time step including control steps
    // to enable it, both extObj and extUpdate have to be set.
	// it can be used for visualization interface
	void* extObj;  // pointer to the external object worked by the hooked function
	void (*extUpdate)(void*, ScadaGenerator*); //function pointer to visual update function

};


enum SG_ERR {  /* Error codes*/
	SG_OK,

	/* below are errors specific to this (sub)class */
	SG_NET_FILE_ERR,
	SG_NET_PROBLEM,
	SG_NET_INIT_ERR,
	SG_ORIG_NET_WONT_RUN,
	SG_CANT_EST_DB,
	SG_MORE_THAN_ONE_PROVIDER,

	SG_SCADA_INIT_ERR,
	SG_UNKNOWN_CHANNEL,
	SG_HYD_PROBLEM,
	SG_CANT_FILL_SCADA,

	SG_FILE_PATH_TOO_LONG,
	SG_FILE_PATH_NOT_ANSI,
	SG_INVALID_VARIMA,
	SG_VARIMA_N_NOT_MATCH,

	SG_MAPPER_PTR_INVALID,
	SG_MEM_NOT_ALLOCED,
    SG_INIT_DATA_NOT_ENOUGH,

	SG_DUMMY_LAST	};

/// Factory method
/** Load EPANET network inp file. Initialize hydraulic engine.
Initialize Provider (database). Initialize Random number generator */
SG_ERR SG_new(TCHAR* inpfile_path, char const* ini_path, int seed, double cov, ScadaGenerator** sg);
SG_ERR SG_new(TCHAR* inpfile_path, char const* ini_path, int seed, Varima* md, ScadaGenerator** sg);

/// Report error
void SG_ewi(SG_ERR err);

/// Run hydraulic simulation and write to database  
SG_ERR SG_make(ScadaGenerator *sg, Tstamp* t0, int timespan);

/// struct for makeT arguments
struct SG_MakeTArg {
	ScadaGenerator* sg;
	Tstamp* t;
	int timespan;
};

/// Run SG_make - entry function for multi-threading
inline unsigned WINAPI SG_makeT(void* args) {
    SG_MakeTArg* arguments = (SG_MakeTArg*) args;
	return (unsigned)SG_make(arguments->sg, arguments->t, arguments->timespan);
};

/// destroy 
SG_ERR SG_delete(ScadaGenerator *sg);

