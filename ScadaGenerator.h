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
#include <process.h>
#include "DataSource.h"
#include "Varima.h"
#include "vtkPolyDataMapper.h"
#include "vtkCylinderSource.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include <tchar.h>

class ScadaGenerator: public Provider {

	/* A ScadaGenerator fills the database with generated SCADA data 
from a network file and pre-defined temporal and/or spatial correlations
of water demands */

public:  /* public definitions */
	enum Err {  /* Error codes*/
		OK,
		ENV_NOT_ALLOC,  //can not allocate environment handle
		ENV_CANT_SET_V3,  //can not set env handle to v3
		DBC_NOT_ALLOC,  //can not allocate db connection handle
		CANT_CONNECT_DB,
		STMT_NOT_ALLOC,  //can not allocate statement handle
		CANT_LOAD_CHANNELS,
		UNKNOWN_MSMT_TYPE,
		DB_CONN_NOT_READY,
		CANT_QUERY_DATA,
		NO_DATA_FOR_THE_CHANNEL,
		CANT_PREPARE_QUERY,
		DSN_TOO_LONG_OR_NULL,
		
		/* below are errors specific to this (sub)class */
		NET_FILE_ERR = 0x100,
		NET_PROBLEM,
		ORIG_NET_WONT_RUN,
		UNKNOWN_CHANNEL,
		CHAN_TABLE_ERR,
		HYD_PROBLEM,
		CANT_PREPARE_SQL,
		CANT_FILL_SCADA,
		FILE_PATH_TOO_LONG,
		FILE_PATH_NOT_ANSI,
		INVALID_VARIMA,
		VARIMA_N_NOT_MATCH,
		MAPPER_PTR_INVALID
		
			};

	enum Mode { /* Data generation modes*/
		// Possible randomization scenarios:
		// independent log-normal rv for each node at each step
		INDEPENDENT,
		// temporally and spatially correlated
		TEMPORAL_SPATIAL_CORRELATED
	};

	struct Makeargs {
		ScadaGenerator* sgtr;
		Tstamp* t;
		int timespan;
		Varima* md;
		vtkPolyData** arr; 
		/* which step is data ready */
		volatile LONG* datastep;
		int* dispstep;
	};


public: /* public methods */

	/* init hydraulic engine and check/reset the database tables,
	create sketch of network visualization*/
	Err init(LPCTSTR inpfile_path, 
		LPCTSTR dsn, /* if dsn==NULL don't export to database*/
		int seed);
	Err init(LPCTSTR inpfile_path, 
		LPCTSTR dsn, /* if dsn==NULL don't export to database*/
		int seed, 
		vtkPolyData* nettop[2] /* if mapper==NULL no visualization*/
		);

	/* get number of water users in the network*/
	int getNuser() {return nuser;};
	int getNTstep() {return ntstep;}; // total steps including control steps
	int getDur() {return dur;}; //network sim duration time

	/* generate Scada data for a time range*/
	/* overload function 1, mode INDEPENDENT */
	Err make(
		SQL_TIMESTAMP_STRUCT* t1,  /* The time struct for starting time */
		int timespan, /*The total time of data generation, in seconds*/ 
		double cov /* coefficient of variation for the log-normal distribution */
		);

	/* overload function 2, mode TEMPORAL_CORRELATED */
	Err make(SQL_TIMESTAMP_STRUCT* t1, int timespan, Varima* md);

	Err make(Tstamp* t1, int timespan, Varima* md, 
		/* array of computed polydatas for pressure and demands visualization,
		must have 2*NTsteps poly allocated*/
		vtkPolyData** arr, 
		/* which step is data ready */
		volatile LONG* datstep,
		/* display time in sec */
		int* dispstep);

	//The following two functions make the generation in a seperate thread
	static unsigned WINAPI make(void* args_in) {
		Makeargs* args = (Makeargs*)args_in;
		return (unsigned)static_cast<ScadaGenerator*>(args->sgtr)->make(
			args->t, args->timespan, args->md, 
			args->arr, args->datastep, args->dispstep);
	};  // called by _beginthreadex

		/* run make in a new thread */
	Err makeTd(Tstamp* t1, int timespan, Varima* md, 
		vtkPolyData** arr, 	volatile LONG* datstep, int* dispstep) {
		Makeargs* args = (Makeargs*)malloc(1*sizeof(Makeargs));
		args->sgtr = this;
		args->t = t1; 
		args->timespan = timespan;
		args->md = md;
		args->arr = arr;
		args->datastep =  datstep;
		args->dispstep = dispstep;
		_beginthreadex(NULL, 0, &ScadaGenerator::make, args, 0, NULL);
		return OK;
	}

	ScadaGenerator();
	~ScadaGenerator();

private:

	Mode mode;  /* generation mode */

	int write_db; /*if write to a database */

	long dur; /* epanet simulation duration in seconds */
	long hstep; /* epaent hydraulic time step length */
	int nnode, ntank, njunc; /* counts of nodes, tanks, junctions*/
	int nuser;  //no. of water users, water user must have demand >0 
				//at at least one time step
	int ntstep; //total time steps/snapshots, including control steps
	double minx, maxx, miny, maxy, minz, maxz, mind, maxd; 
	//bands of variables to be visualized


	/* max number of characters in a file path*/
	static const int MAX_FILE_PATH = 255; 

	Channel* chanlist;
	unsigned n_chan;

	/* SQL strings */
	//TCHAR sql_clean[MAX_SQL_LEN];
	//TCHAR sql_create[MAX_SQL_LEN];

	Rvgs* rg;  // random num generator

	double** demands;  //storage for mean demands at all hyd steps
	//demand[i][istep] is the ith nodal demand at timestep istep

	Err make(SQL_TIMESTAMP_STRUCT* t1, int timespan); //make func

	// for independent mode only
	double _cov;

	// for correlated demand mode only
	Varima* _md;
	double** _user;  //water use matrix, stores pointers to water users'
	//  water demands
	double* _tpd; // temporary water demand storage

	// for data visualization only
	int _needsVis;  // need visualization
	vtkPolyData** _net_array;
	volatile LONG* _datastep;
	int* _dispstep;
	vtkTransformFilter *cyl;
	vtkPolyData *zlineP, *zlineD;
	vtkPolyData* _netpoly[2];// reference 3d [0] and 2d [1] networks
	
};




	