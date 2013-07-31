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
#include <Windows.h>
#include <sql.h>
#include <sqlext.h>
#include "Network.h"
#include <tchar.h>
#include <atltime.h>
#include <strsafe.h>

class Provider;
struct Channel;

#define MAX_NET_ID_LEN  255
#define Tstamp SQL_TIMESTAMP_STRUCT  //Use Timestamp struct

class DataSource
/* Data storage and channel/provider manager for the aligned sensor data */
{

public:
	enum Err {OK,
		CANT_ADD_PROVIDER,
		NO_CHAN_ADDED,
		MEM_NOT_ALLOC,
		CANT_GET_DATA_THRU_CHAN
	};
	enum UnitsType {};  // work needed here

	// initialization
	DataSource(unsigned dt_in  /* time quantum in seconds*/);
	DataSource(unsigned dt_in, Provider* prov /* a provider */);

	/* add a data provider and all its channels */
	Err addProvider(Provider*);

	// data retrieval

	/*get number of channels*/
	unsigned getNChan() {return n_chan;};

	/*get a channel with given channel id, if not found return NULL*/
	Channel* getChan(int id);

	// return channel list
	Channel* getChanList() {return channel_list;};

	/*get the n snapshots at or before a specified time, 
	snapshot should already be allocated to contain n_chan*(n+1) doubles*/
	Err getSnapshots(Tstamp t_in, unsigned n, double* snapshots);

	/* get a data snapshot (wrapper) */
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

	/*dump data about channels and providers*/
	void dumpChannelsInfo();


	DataSource();
	~DataSource(); //clean all providers, release memory of the data buffer

protected:
	//static const int BUF_LEN = 1000; // number of time steps record in data buffer

	Channel* channel_list; // list of data channels
	int n_prov, n_chan;  //number of data providers and channels
	unsigned	dt;  // time interval between snapshots, in seconds

};



class Provider
/*  ODBC-interfaced Data provider, stores a database access string. 
	A provider may have one or multiple channels */
{
protected:
	static const int MAX_DSN_LEN = 2048;
	static const int MAX_TABLE_NAME = 1024;
	static const int MAX_SQL_LEN = 1023;
	SQLHENV     hEnv;   //environment handle
    SQLHDBC     hDbc;	//connection handle
    SQLHSTMT    hStmt;	//general statement handle
	SQLHSTMT	hStmt_dat; //statement for getting data
	int		hStmt_dat_prepared;  // flag the hStmt_dat is prepared

public:
	Provider();

	enum Err {
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
		DSN_TOO_LONG_OR_NULL
		
			};
	Err connect(LPCTSTR dsn);  /*odbc data source name*/

	Err check(LPCTSTR dat_table_name,
			 LPCTSTR meta_table_name); /* check data tables ok */ 

	/* create channels and load channel information */
	Err loadChannels(Channel** /* output the list head of channels */
		, unsigned* /*output the number of channels*/); 

	/* get a (observational) data point from a channel 
	at a specific timestamp or left-closest to the timestamp*/
	Err getDataAt(SQL_TIMESTAMP_STRUCT timestamp, 
			  int timeshift,  /* desired time shift of time stamp, in sec*/
			  int chan_key,   /* the desired channel*/
			  double* data_out,    /* The data */
			  int* time_offset_out  /*diff of desired and returned timestamp, 
								in seconds, NULL if this info is not needed*/ 
			  );
	

	/* destructor, disconnect db connection, free channels*/
	~Provider();

protected:
	// interpret odbc errors
	void handleDiagnosticRecord(SQLHANDLE      hHandle,    
                             SQLSMALLINT    hType,  
                             RETCODE        RetCode);

protected:
	static const TCHAR dat_tab[MAX_TABLE_NAME];  // fact table name string
	static const TCHAR chan_tab[MAX_TABLE_NAME];  // metadata table name string
};


struct Channel
	/* a channel is a means through which data collected 
	by a sensor can be accessed.*/
{
	int			key;		//channel key (e.g., database primary key)
	char		name[MAX_NET_ID_LEN];  // name of the component, consistent with network
						// component id string
	Provider*	provider;  //data provider
	Network::FieldType  mtype;  // the type of measurement in a hydraulic model
	int			mindex;  //index in the hydraulic network
	DataSource::UnitsType  unit; //unit of measurement
	double		rse;  // relative standard error,=std. dev/mean,(assuming normal)
	double		lower_lim;  //lower limits of the possible measurements
	double		upper_lim;  //upper limits
	double		err_data;   // if the the sensor has anormaly, this value shows up.

	enum Status {OK, DATA_ERR, DATA_OUTBOUND};
	Status		status;  // working status of this channel

	Channel*	next;   //next channel in a list
};


