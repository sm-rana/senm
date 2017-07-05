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

/* Classes and structs associated with the Database Access Layer*/
#pragma once
#include <Windows.h>
#include <sql.h>
#include <sqlext.h>
#include <tchar.h>
#include <atltime.h>
#include <strsafe.h>
#include "SenmCoreIncs.h"
#include "Channel.h"
#include "Network.h"

//struct Channel;
struct Provider;
//struct DataSource;

#define Tstamp SQL_TIMESTAMP_STRUCT  ///> Use Timestamp struct
#define DATASOURCE_MAX_CHAN_INFO  128  ///> max string size of channel type info text

/// Data source manager
/** The instance of DataSource class manages a list of data channels and provides 
data services. must be created with DataSource::New()
*/
struct DataSource {
	/** Status of data source*/
	enum Err {OK,
		MEM_NOT_ALLOCED,
		PROVIDER_INVALID,
		CANT_ADD_PROVIDER,
		NO_CHAN_ADDED,
		CANT_GET_DATA_THRU_CHAN,
		INI_FILE_ERR,
		DELTA_T_ERR_INI_FILE,
		N_PROVIDER_ERR_INI_FILE,
		CANNOT_CREATE_PROVIDER,

		DUMMY_LAST

	};

	// work needed here, for now assume network and channels always have the same unit 
	enum UnitsType {};  

	int n_prov;  // number of data providers 
	int n_chan;  // number of data channels
	int	dt;  // time interval between snapshots, in seconds

	double* _datBuf;  // data buffer - most recent snapshot with updateBufSnapshot()
	Provider* lsProv; // list of data providers, mem managed here
	Channel* lsChan; // list of channels, mem managed by provider

	/// report text representation of error 
	static void reportEWI(Err err);

	/// Factory method: create a DataSource with an existing data provider. 
	/** \param [in] dt_in     Time quantum in seconds
	\param [in] prov      Initial data provider
	\param [out] ds_out    newly creately datasource
	*/
	//static Err New(int dt_in, Provider* prov, DataSource** ds_out);

	/// Factory method: create a DataSource using a .ini configuration file
	/** \param [in] ininame    INI configuration file name
	\param [out] ds_out   created datasource object
	*/
	static Err New(const char *ininame, DataSource** ds_out);
	static Err New(const char *ininame, Network* net, DataSource** ds_out);

	///
	/** add a data provider and all its channels */
	Err addProvider(Provider*);

	/// update databuffer with the snapshot at a given time.
	Err updateBufSnapshot(CTime ct);

	/** Print information about channels and providers*/
	void dumpChannelsInfo();




	/** Get a channel with given channel id, if not found return NULL*/
	/** \return     a pointer to a channel, NULL if it does not exist*/
	Channel * getChan(int id);

	/// Get consecutive n snapshots from the channels into pre-allocated arrays.
	/**The method gets n consecutive snapshots at and before a specified time. 
	Snapshots are taken at t_in, t_in-dt, t_in-2*dt, ...t_in - (n_ss-1)*dt
	The pointer snapshots should have n_chan*(n_ss) doubles allocated already.
	\param [in]  t_in    Time stamp on which the data query is based	
	\param [in]  n_ss       Number of snapshots needed 
	\param [out]  snapshots   pointer to a double array
	*/
	Err fillSnapshots(Tstamp t_in, int n_ss, double* snapshots);

	/// Get a single snapshot of data into pre-allocated arrays.
	/** Get the most recent snapshot as of ct
	\param [in] ct      CTime struct representing time in need.
	\param [out] snapshot   Pointer to an array, should have been allocated with n_chan doubles
	*/
	Err fillASnapshot(CTime ct, double* snapshot);

	~DataSource(); //clean all providers, release memory of the data buffer

protected: //disable default constructor
	DataSource();


};



#define MAX_DSN_LEN  2048 ///> Maximum length of the DSN (connection string)
#define MAX_TABLE_NAME  1024 ///> Maximum size of table name string
#define MAX_SQL_LEN  1023 ///> Maximum size of SQL query string

/// Database connection manager 
/**  An instance of the Provider class manages a database connection.
The db query funtions are implemented using ODBC.
-interfaced Data provider, stores a database access string. 
A provider may have one or multiple channels */
struct Provider
{
	// TODO: support multiple database implementations
	SQLHENV     hEnv;   ///> environment handle
	SQLHDBC     hDbc;	///> connection handle
	SQLHSTMT    hStmt;	///> general statement handle
	SQLHSTMT	hStmt_dat; ///> statement for reading data
	SQLHSTMT	hStmt_w; ///> statement for writing data (scada generator)
	int		hStmt_dat_prepared;  ///>  flag showing if the hStmt_dat handle is ready (prepared)
	Channel*  lsChan; ///> channel list 
	int nChan; ///> number of channels

	/** Data connection error state*/
	enum Err {
		OK,
		MEM_NOT_ALLOCED, 

		ENV_NOT_ALLOC,  ///>can not allocate environment handle
		ENV_CANT_SET_V3,  ///>can not set env handle to v3
		DBC_NOT_ALLOC,  ///>can not allocate db connection handle
		DSN_TOO_LONG_OR_NULL,
		CANT_CONNECT_DB,
		STMT_NOT_ALLOC,  ///>can not allocate statement handle

		CANT_PREPARE_SQL,
		CANT_INSERT_SCADA,
		CANT_LOAD_CHANNELS,
		UNKNOWN_MSMT_TYPE,
		DB_CONN_NOT_READY,

		CANT_QUERY_DATA,
		NO_DATA_FOR_THE_CHANNEL,
		CANT_PREPARE_QUERY,
		CANT_RESET_FACT_TAB,

        CANT_FIND_COR_COMP,

		DUMMY_LAST
	};

	/// Error reporting
	static void reportEWI(Err err);

	/// Factory method
	/** Connect to the given database and load channel metadata
	   if an optional Network is provided, mindex field of the channels
	   will have correct values*/
	static Err New(TCHAR* tdsn, Provider** prov_out, Network* net = NULL);

	/**  Connect to a database */
	/** There must be two data tables existing in the db.
	- Msmts  Measurement data (fact table)
	- Channels   Channel information used to construct the list of channels
	*/
	/** \param [in] dsn      odbc database connection string
	*/
	Err connect(LPCTSTR tdsn);  

	// TODO: Diagnostic checking of data tables */
	//Err check(LPCTSTR dat_table_name,
	//		 LPCTSTR meta_table_name); 

	/** Create channels and load channel information */
	/**  Load channel information from the "Channels" table in the database
        if net is not null, mindex of the channels will be looked up
	*/
	Err loadChannels(Network* net); 

	/// Insert Scada records - initialization
	Err initInsertScada(Tstamp* scada_t0);

	/** Insert scada records*/
	Err insertScada(int* timeshift, int* chan_key, float* signal);

	/** Reset fact table, clean all scada data*/
	Err resetFactTab();

	/** Get an observed data point from a channel by db primary key 
	/**  Obtain a data point at a specified time, if there is no matching time,
	return the closest data point observed right before that time.
	\param [in] timestamp    Timestamp
	\param [in] timeshift    Shift back this amount of time to the timestamp
	\param [in] chan_key     Primay Key (cid) in the Channel table, identify a channel
	\param [out] data_out    pointer to which the data retrieved will be stored
	\param [out] time_offset_out   pointer where the function puts the time in seconds between 
	the requested timestamp and the return timestamp (in case there is no matching time).
	if this pointer == NULL, nothing will be returned.
	*/
	Err getDataAt(Tstamp timestamp, 
		int timeshift,  /* desired time shift of time stamp, in sec*/
		int chan_key,   /* the desired channel*/
		double* data_out,    /* The data */
		int* time_offset_out  /*diff of desired and returned timestamp, 
							  in seconds, NULL if this info is not needed*/ 
							  );


	/** destructor, 
	/**disconnect db session*/
	~Provider();

	/// interpret odbc errors
	void handleDiagnosticRecord(SQLHANDLE      hHandle,    
		SQLSMALLINT    hType,  
		RETCODE        RetCode);

	TCHAR dat_tab[MAX_TABLE_NAME];  ///> fact table name string, set to "Msmts"
	TCHAR chan_tab[MAX_TABLE_NAME];  ///> metadata table name string, set to "Channels"

	Provider* next ; ///> next provider in a list

protected:
	Provider();
};
