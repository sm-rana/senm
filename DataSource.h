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
#include "Network.h"
#include <tchar.h>
#include <atltime.h>
#include <strsafe.h>

class Provider;
struct Channel;

#define MAX_NET_ID_LEN  255 ///> Network node/link identification string
#define Tstamp SQL_TIMESTAMP_STRUCT  ///> Use Timestamp struct

/// Data source manager
/** The instance of DataSource class manages a list of data channels and provides 
   data services.
   */
class DataSource
{

public:
    /** Status of data source*/
	enum Err {OK,
		CANT_ADD_PROVIDER,
		NO_CHAN_ADDED,
		MEM_NOT_ALLOC,
		CANT_GET_DATA_THRU_CHAN
	};
	enum UnitsType {};  // work needed here

	/// Constructor. Data connection needs to be added later via addProvider(). 
    /// \param [in] dt_in     Time quantum in seconds
	DataSource(unsigned dt_in );

    /// Construct the data source with an existing data provider. 
    /** \param [in] dt_in     Time quantum in seconds
        \param [in] prov      Initial data provider
        */
	DataSource(unsigned dt_in, Provider* prov /* a provider */);

	/** add a data provider and all its channels */
	Err addProvider(Provider*);

	// @{

	/**get number of channels*/
	unsigned getNChan() {return n_chan;};

	/** Get a channel with given channel id, if not found return NULL*/
    /** \return     a pointer to a channel, NULL if it does not exist*/
	Channel* getChan(int id);

	/** Get the list of channels */
	/** \return    channel list */
	Channel* getChanList() {return channel_list;};

    /** Get consecutive n snapshots from the channels.
	The method gets n consecutive snapshots at and before a specified time. 
    Snapshots are taken at t_in, t_in-dt, t_in-2*dt, ...
	The pointer snapshots should have n_chan*(n+1) doubles allocated already.
      \param [in]  t_in    Time stamp on which the data query is based	
      \param [in]  n       Number of snapshots needed - 1
      \param [out]  snapshots   pointer to a double array
	*/
	Err getSnapshots(Tstamp t_in, unsigned n, double* snapshots);

	/** Get a single data snapshot (wrapper) */
    /** Get the most recent snapshot as of ct
      \param [in] ct      CTime struct representing time in need.
      \param [out] snapshot   Pointer to an array, should have been allocated with n_chan doubles
      */
	Err getASnapshot(CTime ct, double* snapshot){
		Tstamp tp;
		tp.year = ct.GetYear();
		tp.month = ct.GetMonth();
		tp.day = ct.GetDay();
		tp.hour = ct.GetHour();
		tp.minute = ct.GetMinute();
		tp.second = ct.GetSecond();
		tp.fraction = 0;
		return getSnapshots(tp, 0, snapshot);
	};

	/** Print information about channels and providers*/
	void dumpChannelsInfo();

    // @}


	DataSource();
	~DataSource(); //clean all providers, release memory of the data buffer

protected:
	//static const int BUF_LEN = 1000; // number of time steps record in data buffer

	Channel* channel_list; ///> list of data channels
	int n_prov;  ///> number of data providers 
    int n_chan;  ///> number of data channels
	unsigned	dt;  ///> time interval between snapshots, in seconds

};


/// Database connection manager 
/**  An instance of the Provider class manages a database connection.
    The db query funtions are implemented based on ODBC.
	-interfaced Data provider, stores a database access string. 
	A provider may have one or multiple channels */
class Provider

{
protected:
	static const int MAX_DSN_LEN = 2048; ///> Maximum length of the DSN (connection string)
	static const int MAX_TABLE_NAME = 1024; ///> Maximum size of table name string
	static const int MAX_SQL_LEN = 1023; ///> Maximum size of SQL query string
	SQLHENV     hEnv;   ///> environment handle
    SQLHDBC     hDbc;	///> connection handle
    SQLHSTMT    hStmt;	///> general statement handle
	SQLHSTMT	hStmt_dat; ///> statement for getting data
	int		hStmt_dat_prepared;  ///>  flag showing if the hStmt_dat handle is ready (prepared)

public:
	Provider();
    
    /** Data connection error state*/
	enum Err {
		OK,
		ENV_NOT_ALLOC,  ///>can not allocate environment handle
		ENV_CANT_SET_V3,  ///>can not set env handle to v3
		DBC_NOT_ALLOC,  ///>can not allocate db connection handle
		CANT_CONNECT_DB,
		STMT_NOT_ALLOC,  ///>can not allocate statement handle
		CANT_LOAD_CHANNELS,
		UNKNOWN_MSMT_TYPE,
		DB_CONN_NOT_READY,
		CANT_QUERY_DATA,
		NO_DATA_FOR_THE_CHANNEL,
		CANT_PREPARE_QUERY,
		DSN_TOO_LONG_OR_NULL
		
			};

    /**  Connect a database */
    /** There must be two data tables exsited in the db.
       - Msmts  Measurement data (fact table)
       - Channels   Channel information used to construct the list of channels
       */
    /** \param [in] dsn      odbc database connection string
    */
	Err connect(LPCTSTR dsn);  

    /** Diagnose check data tables */
    /** Not implemented yet*/
	Err check(LPCTSTR dat_table_name,
			 LPCTSTR meta_table_name); 

	/** Create channels and load channel information */
    /**  Load channel information from the "Channels" table in the database
         \param [out] chan_list     a pointer where the pointer to the Channel list will be output
         \param [out] n_chan        a pointer where the number of channel retrieved be output 
         */
	Err loadChannels(Channel** chan_list/* output the list head of channels */
		, unsigned* n_chan/*output the number of channels*/); 

	/** Get an observed data point from a channel 
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

protected:
	/// interpret odbc errors
	void handleDiagnosticRecord(SQLHANDLE      hHandle,    
                             SQLSMALLINT    hType,  
                             RETCODE        RetCode);

protected:
	static const TCHAR dat_tab[MAX_TABLE_NAME];  ///> fact table name string, set to "Msmts"
	static const TCHAR chan_tab[MAX_TABLE_NAME];  ///> metadata table name string, set to "Channels"
};

/// Data channel
 /** The struct stores descriptive information about a data channel.
 Data channels are created by Provider::loadChannels(),
 and deleted by DataSource::~Datasource()
 */
struct Channel
{
	int			key;		///> channel key (e.g., database primary key)
	char		name[MAX_NET_ID_LEN];  ///>  name of the component, consistent with network component id string
	Provider*	provider;  ///> data provider
	Network::FieldType  mtype;  ///>  the type of measurement in a hydraulic model
	int			mindex;  ///> index in the hydraulic network
	DataSource::UnitsType  unit; ///> unit of measurement
	double		rse;  ///>  relative standard error,=std. dev/mean,(assuming normal)
	double		lower_lim;  ///> lower limits of the possible measurements
	double		upper_lim;  ///> upper limits
	double		err_data;   ///>  if the the sensor has anormaly, this value shows up.

	enum Status {OK, DATA_ERR, DATA_OUTBOUND};
	Status		status;  ///>  working status of this channel

	Channel*	next;   ///> next channel in a list
};


