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
#include "DataSource.h"

class ScadaGenerator: public Provider {

	/* A ScadaGenerator fills the database with generated SCADA data 
from a network file and pre-defined temporal and/or spatial correlations
of water demands */

public:
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
		
		/* below are errors specific to the subclass */
		NET_FILE_ERR = 1024,
		NET_PROBLEM,
		ORIG_NET_WONT_RUN,
		UNKNOWN_CHANNEL,
		CHAN_TABLE_ERR,
		HYD_PROBLEM,
		CANT_PREPARE_SQL,
		CANT_FILL_SCADA
		
			};

	enum Mode { /* Data generation modes*/
		// Possible randomization scenarios:
		// independent log-normal rv for each node at each step
		INDEPENDENT,
		// temporally correlated demands (e.g, log-AR(2))
		TEMPORAL_CORRELATED,
		// spatially correlated demands (
		SPATIAL_CORRELATED,
		// temporally and spatially correlated
		TEMPORAL_SPATIAL_CORRELATED
	};


	/* empty constructor */
	ScadaGenerator();

	/* init hydraulic engine and check/reset the database tables*/
	Err init(const char* inpfile_path, const char *dsn);

	/* generate Scada data for a time range*/
	/* overload function 1, mode INDEPENDENT */
	Err make(SQL_TIMESTAMP_STRUCT* t1, int timespan /*in seconds*/, double cof);
	
	/* overload function 2, mode TEMPORAL_CORRELATED */
	//Err make(SQL_TIMESTAMP_STRUCT* t1, int timespan, ARIMAmodel md);

	~ScadaGenerator();

private:

	long dur; /* epanet simulation duration in seconds */
	long hstep; /* epaent hydraulic time step length */
	int nnode, ntank, njunc; /* counts of nodes, tanks, junctions, and controls*/

	static const int MAX_FILE_PATH = 255;


	Channel* chanlist;
	unsigned n_chan;

	/* SQL strings */
	
	char sql_clean[MAX_SQL_LEN];
	char sql_create[MAX_SQL_LEN];
	double** demands;  //storage for mean demands at all hyd steps

};




	