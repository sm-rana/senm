#include <Windows.h>
#include <sqlext.h>
#include <atlconv.h>
#include "iniparser.h"
#include "SenmCoreIncs.h"
#include "DataSource.h"

#ifdef _DEBUG
#define debug 1
#else
#define debug 0
#endif

/* Class DataSource implementation */
DataSource::DataSource() {
	_datBuf = NULL;
	lsChan = NULL;
	lsProv = NULL;
	n_prov = 0; 
	n_chan = 0;
	dt = 0;
}

DataSource::Err DataSource::New(const char * ininame, DataSource** ds) {
    return New(ininame, NULL, ds);
}

DataSource::Err DataSource::New(const char * ininame, Network* net, DataSource** ds) {
	USES_CONVERSION;  // for using A2T macro
	*ds = new DataSource();
	if (*ds== NULL) return DataSource::MEM_NOT_ALLOCED;

	Err err;
	dictionary* dict;
	dict = iniparser_load(ininame);
	if (dict == NULL) {
		err = INI_FILE_ERR;
		goto CLEAN; 
	}

	int dt = iniparser_getint(dict, "General:deltat", -1);
	if (dt <= 0 ) {
		err = DELTA_T_ERR_INI_FILE;
		goto CLEAN;
	} else (*ds)->dt = dt;

	int nprov = iniparser_getint(dict, "General:providers", 0);
	if (nprov <= 0) {
		err = N_PROVIDER_ERR_INI_FILE;
		goto CLEAN;
	}; //else (*ds)->n_prov = nprov;

	for (int iprov=1; iprov<=nprov; ++iprov) {
		char tmp[MAX_DSN_LEN];
		char *dsn, *server, *port, *uid, *pwd ;
		char *vendor, *schema;

		sprintf(tmp, "Source%d:vendor", iprov);
		vendor = iniparser_getstring(dict, tmp, "");

		sprintf(tmp, "Source%d:dsn", iprov);
		dsn = iniparser_getstring(dict, tmp, "");

		sprintf(tmp, "Source%d:server", iprov);
		server = iniparser_getstring(dict, tmp, "");

		sprintf(tmp, "Source%d:port", iprov);
		port = iniparser_getstring(dict, tmp, "");

		sprintf(tmp, "Source%d:uid", iprov);
		uid = iniparser_getstring(dict, tmp, "");

		sprintf(tmp, "Source%d:pwd", iprov);
		pwd = iniparser_getstring(dict, tmp, "");

		char src[MAX_DSN_LEN];
		sprintf(src, "DSN=%s; SERVER=%s; PORT=%s; UID=%s; PWD=%s", 
			dsn, server, port, uid, pwd);

		if (strcmp(vendor, "postgres") == 0) {
			sprintf(tmp, "Source%d:schema", iprov);
			schema = iniparser_getstring(dict, tmp, "");
			strcat(src, "; A6=set search_path to ");
			strcat(src, schema);
			strcat(src, ",public;");
		}

		TCHAR tdsn[MAX_DSN_LEN];
		_tcscpy(tdsn, A2T(src));

		Provider* aPrv;
		if (Provider::Err perr = Provider::New(tdsn, &aPrv, net)) {
			Provider::reportEWI(perr);
			err = CANNOT_CREATE_PROVIDER;
			goto CLEAN;
		}

		if (Err derr = (*ds)->addProvider(aPrv)) {
			reportEWI(derr);
			err = derr;
			goto CLEAN;
		}
	}
	return OK;

CLEAN:
	iniparser_freedict(dict);
	delete (*ds);
	(*ds) = NULL;
	return err;

}

void DataSource::reportEWI(Err err) {
	TCHAR errtxt[DUMMY_LAST][MAX_ERR_STRING_SIZE] = {
		TEXT("OK."),
		TEXT("Could not allocate memory for an object or array."),
		TEXT("Data provider provided is not valid."),
		TEXT("Could not add a provider."),
		TEXT("No channel is added."),
		TEXT("Could not obtain data from the channels."),
		TEXT("Could not load INI config file."),
		TEXT("Error reading delta t in ini file."),
		TEXT("Error reading number of providers in ini file."),
		TEXT("Could not create provider using the given config file"),

	};

	TCHAR txt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];
    _stprintf(txt, TEXT("DataSource:%s"), errtxt[err]);

    ewi(txt);
}



DataSource::Err DataSource::addProvider(Provider* prov_in) {
	unsigned n_new_chan;
	//Provider::Err perr;

	if (prov_in == NULL) return PROVIDER_INVALID;
//	if (perr = prov_in->loadChannels(NULL)) {
//		Provider::reportEWI(perr);
//		return CANT_ADD_PROVIDER;
//	}
	n_new_chan = prov_in->nChan;
	if ( n_new_chan == 0) {
		return NO_CHAN_ADDED;
	}

	Channel* iChan = prov_in->lsChan;
	while (iChan->next) iChan = iChan->next; // go to tail of the new channel list to be added
	iChan->next = lsChan; // 
	lsChan = prov_in->lsChan; //update channel list head

	// increase buffer size
	double* _tpdatBuf = (double*) calloc(n_chan + n_new_chan, sizeof(double));
	if (_tpdatBuf == NULL) return MEM_NOT_ALLOCED;
	for (int iich=0; iich<n_chan; ++iich) _tpdatBuf[iich] = _datBuf[iich];

	free(_datBuf);

	// insert the new provider
	prov_in->next = lsProv;
	lsProv = prov_in; // update provider list head

	n_chan += n_new_chan;
	_datBuf = _tpdatBuf;
	++ n_prov;
	return OK;
}

Channel* DataSource::getChan(int chan_id) {
	for (Channel* chan_it=lsChan; chan_it; chan_it = chan_it->next) 
		if (chan_it->key == chan_id) return chan_it;
	return NULL;
}

void DataSource::dumpChannelsInfo() {
	USES_CONVERSION;
	TCHAR txts[][DATASOURCE_MAX_CHAN_INFO] = {
		TEXT("No channel"),
		TEXT("Water level"),
		TEXT("Pump flow rate"),
		TEXT("Minor headloos coefficient"),
		TEXT("Discharge-side pressure"),
		TEXT("Suction-side pressure"),
		TEXT("Pipe flow rate "),
		TEXT("Pipe/valve on/off"),
		TEXT("Nodal pressure"),
		TEXT("Real-time demad (flow meter)"),
        TEXT("Unknown measurement")
	};

	_ftprintf(stderr, 
		TEXT("DataSource Info\nTime Quantum: %d Sec, Total %d channels from %d providers.\n"), 
		dt, n_chan, n_prov);

	TCHAR* type, *status;

	for(Channel* chan_it=lsChan; chan_it; chan_it = chan_it->next) {
		switch (chan_it->type) {
		case Channel::L: type = txts[1]; break;
		case Channel::F: type = txts[2]; break;
		case Channel::V: type = txts[3]; break;
		case Channel::B: type = txts[4]; break;
		case Channel::A: type = txts[5]; break;
		case Channel::Q: type = txts[6]; break;
		case Channel::C: type = txts[7]; break;
		case Channel::P: type = txts[8]; break;
		case Channel::D: type = txts[9]; break;
		default: type= txts[10];
		}

		switch (chan_it->status) {
		case Channel::OK: status = TEXT("OK"); break;
		case Channel::DATA_ERR: case Channel::DATA_OUTBOUND:
			status = TEXT("Anormaly detected"); break;
		}

		_ftprintf(stderr, TEXT(
			"Channel {%d}: %s @ EPANET Idx %d(%s): stderr:%g\n"), 
			chan_it->key, type, chan_it->mindex, A2T(chan_it->name), chan_it->stde);
	}
}

DataSource::Err DataSource::fillSnapshots(Tstamp t_in, int n, double* snapshots) {
    Err err = OK;
	//snapshots memory must be alloced,
	if (snapshots == NULL) return MEM_NOT_ALLOCED;
	int j = 0;  //j - channel no.
	int fErr = 0;

	for (int i=0; i<n; ++i) { // iterate through time steps
		j=0;
		for (Channel* chan_it=lsChan; chan_it; chan_it=chan_it->next) { 
				Provider::Err err = chan_it->provider->getDataAt(
					t_in, -(n-1-i)*dt, /* time shifts backward*/
					chan_it->key, &snapshots[i*n_chan + (j++)], NULL);
				if (err) {
					Provider::reportEWI(err);
					fErr = 1;
				}
		}
	}

	if (fErr) {
        reportEWI(err = CANT_GET_DATA_THRU_CHAN);
	} 
    return err;
}

DataSource::Err DataSource::fillASnapshot(CTime ct, double* snapshot) {
	Tstamp tp;
	tp.year = ct.GetYear();
	tp.month = ct.GetMonth();
	tp.day = ct.GetDay();
	tp.hour = ct.GetHour();
	tp.minute = ct.GetMinute();
	tp.second = ct.GetSecond();
	tp.fraction = 0;
	Err ec =  fillSnapshots(tp, 0, snapshot);
	//update data buffer
	//    for (int iChan = 0; iChan<n_chan; ++iChan) {
	//      _datBuf[iChan] = snapshot[iChan];
	//	}
	return ec;
}

DataSource::Err DataSource::updateBufSnapshot(CTime ct) {
	Err ec = fillASnapshot(ct, _datBuf);
	return ec;
}


DataSource::~DataSource() {
	//destroy all channels
	//	Channel* head = channel_list;
	//	Channel* cur;
	//	while (head) {
	//		cur = head;
	//		head = head->next;
	//		delete cur;
	//	}
	free(_datBuf);

	//destroy providers (which in turn destroy channels)

	Provider* tpPv = lsProv;
	Provider* cur;
	while (tpPv) {
		cur = tpPv;
		tpPv = tpPv->next;
		delete cur;
	}
}

/* Class Provider implementation */

Provider::Provider() 
{
	hEnv = NULL;
	hDbc = NULL;
	hStmt = NULL;
	hStmt_dat_prepared = 0;
	lsChan = NULL;
	nChan = 0;
	next = NULL;
	_tcscpy(dat_tab, TEXT("msmts"));  // fact table
	_tcscpy(chan_tab, TEXT("channels"));  // data table

}

Provider::Err Provider::New(TCHAR* tdsn, Provider** prov_out, Network* net) {
	(*prov_out)  = new Provider();
	if (*prov_out == NULL) return MEM_NOT_ALLOCED;

	Err err;
	err = (*prov_out)->connect(tdsn);
	if (err)  goto END;

	err = (*prov_out)->loadChannels(net);
	if (err) goto END;

	return OK;

END:
	delete (*prov_out);
	return err;

}

void Provider::reportEWI(Err err) {
	TCHAR errtxt[DUMMY_LAST][MAX_ERR_STRING_SIZE] = {
		TEXT("OK"),
		TEXT("Memory can not be allocated."),

		TEXT("Unable to allocate an environment handle"),
		TEXT("can not set env handle to v3"),
		TEXT("can not allocate db connection handle"),
		TEXT("DSN string is either null or too long."),
		TEXT("Cannot connect to database. Check the DSN. "),
		TEXT("Unable to allocate statement handle"),

		TEXT("Could not prepare the SQL statement."),
        TEXT("Could not insert scada data into the database."),
		TEXT("Cannot load data channels. check if table and column exist in the db."),
		TEXT("the type of the measurments is unkonwn."),
		TEXT("Database, statment, or environement handle is not ready. run connect() first."),

		TEXT("SQL statement fails to execute. Do the required columns exist in the table?"),
		TEXT("SQL statement fails to fetch"),
		TEXT("The SQL statment is not prepared. Do the data table/view exist in the db?"),
        TEXT("Could not reset fact table"),
	};

	TCHAR txt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];

	_stprintf(txt, TEXT("<%s>DataProvider:%s"), TEXT("DAL"), errtxt[err]);

    ewi(txt);
}


// Connect to a database with a given dsn string
Provider::Err Provider::connect(LPCTSTR dsn_in  /* Data source name*/ ) {

	SQLTCHAR	dsn[MAX_DSN_LEN];  // data source name
	RETCODE		rc;  //return code

	// Allocate an environment

	if (SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &hEnv) == SQL_ERROR)
	{
		return ENV_NOT_ALLOC;
	}

	// set odbc version
	if (rc = SQLSetEnvAttr(hEnv,
		SQL_ATTR_ODBC_VERSION,
		(SQLPOINTER)SQL_OV_ODBC3,
		0)) {
			//not sucessfuly
			handleDiagnosticRecord(hEnv, SQL_HANDLE_ENV, rc);
			if (rc==SQL_ERROR) return ENV_CANT_SET_V3;
	}

	// allocate connection handle
	if (rc = SQLAllocHandle(SQL_HANDLE_DBC, hEnv, &hDbc)) {
		//if not ok
		handleDiagnosticRecord(hEnv, SQL_HANDLE_ENV, rc);
		if (rc==SQL_ERROR) return DBC_NOT_ALLOC;
	}

	// connect to the database
	// make a copy of dsn
	if (_tcscpy_s(dsn, MAX_DSN_LEN, dsn_in)) // DSN TOO LONG
		return DSN_TOO_LONG_OR_NULL;

	SQLTCHAR dsn_out[MAX_DSN_LEN];

	if (rc = SQLDriverConnect(hDbc, /*Connection*/
		GetDesktopWindow(), /*no window handle is applicable*/
		dsn, /*full connection string*/
		SQL_NTS,  /*null-terminated*/
		dsn_out, MAX_DSN_LEN, NULL, /* don't output dsn*/
		SQL_DRIVER_COMPLETE 
		/* don't prompt user for confidential if info provided in dsn*/
		)) {
			//connection problem
			handleDiagnosticRecord(hDbc, SQL_HANDLE_DBC, rc);
			if (rc==SQL_ERROR) return CANT_CONNECT_DB;
	}// else {
	//	_ftprintf(stderr, TEXT("DB connected.\n"));
	//}

	//Allocate the statement handle
	if (rc = SQLAllocHandle(SQL_HANDLE_STMT, hDbc, &hStmt)) {
		//if not good
		handleDiagnosticRecord(hDbc, SQL_HANDLE_DBC, rc);
		if (rc == SQL_ERROR) return STMT_NOT_ALLOC;
	}
	if (rc = SQLAllocHandle(SQL_HANDLE_STMT, hDbc, &hStmt_dat)) {
		//if not good
		handleDiagnosticRecord(hDbc, SQL_HANDLE_DBC, rc);
		if (rc == SQL_ERROR) return STMT_NOT_ALLOC;
	}
	if (rc = SQLAllocHandle(SQL_HANDLE_STMT, hDbc, &hStmt_w)) {
		//if not good
		handleDiagnosticRecord(hDbc, SQL_HANDLE_DBC, rc);
		if (rc == SQL_ERROR) return STMT_NOT_ALLOC;
	}


	return OK;

}

Provider::Err Provider::initInsertScada(Tstamp* scada_t0){

	SQLTCHAR sqltext[MAX_SQL_LEN];
	RETCODE rc;
    Err err;

	if (debug) _stprintf(sqltext, TEXT(
		"INSERT INTO %s (time, value, cid) VALUES( ?::timestamp + ((? || ' SECONDS')::interval), ?, ?);"),  //postgres
		//"		VALUES(DATE_ADD(?, INTERVAL ? SECOND), ?, ?);"), //mysql
		dat_tab);

	rc = SQLPrepare(hStmt_w, sqltext, SQL_NTS);
	if (rc) {
		handleDiagnosticRecord(hStmt_w, SQL_HANDLE_STMT, rc);
		if (rc == SQL_ERROR) {
			reportEWI(err = CANT_PREPARE_SQL); return err;
		}
	}

	rc = SQLBindParameter(hStmt_w, 1, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, scada_t0, 0, NULL);
    return OK;
}

Provider::Err Provider::insertScada(int* timeshift, int* chan_key, float* signal) {
	RETCODE rc;
    Err err;

	rc = SQLBindParameter(hStmt_w, 2, SQL_PARAM_INPUT, 
		SQL_C_ULONG, SQL_INTEGER, 10, 0, timeshift, 0, NULL);

	rc = SQLBindParameter(hStmt_w, 4, SQL_PARAM_INPUT, 
		SQL_C_ULONG, SQL_INTEGER, 10, 0, (chan_key), 0, NULL);
	rc = SQLBindParameter(hStmt_w, 3, SQL_PARAM_INPUT, 
		SQL_C_FLOAT, SQL_FLOAT, 15, 0, signal, 0, NULL);

	rc = SQLExecute(hStmt_w);
	if (rc) {
		handleDiagnosticRecord(hStmt_w, SQL_HANDLE_STMT, rc);
		if (rc == SQL_ERROR) {
			reportEWI(err = CANT_INSERT_SCADA);
            return err;
		}
	}
    return OK;
}

Provider::Err Provider::resetFactTab(){
	if ((hEnv == NULL) || (hDbc = NULL) || (hStmt == NULL) )
		return DB_CONN_NOT_READY;

	SQLTCHAR sqltext[MAX_SQL_LEN];
	_stprintf(sqltext, TEXT("DELETE from %s"), dat_tab);

    SQLRETURN rc;
	rc = SQLExecDirect(hStmt, sqltext, SQL_NTS);
    if (rc) {
			handleDiagnosticRecord(hStmt_dat, SQL_HANDLE_STMT, rc);
			if (rc==SQL_ERROR) return CANT_RESET_FACT_TAB;
	}

    return OK;
	
}



Provider::Err Provider::getDataAt(SQL_TIMESTAMP_STRUCT timestamp, 
								  int timeshift,  /* desired time shift of time stamp, in sec*/
								  int chan_key,   /* the desired channel*/
								  double* data_out,    /* The data */
								  int* time_offset /*timestamp of the data returned */ 
								  ) 
{
	if ((hEnv == NULL) || (hDbc == NULL) || (hStmt == NULL) )
		return DB_CONN_NOT_READY;

	SQLTCHAR sqltext[MAX_SQL_LEN];
	SQLTCHAR postgres_time_diff_sql[MAX_SQL_LEN];
	SQLTCHAR mysql_time_diff_sql[MAX_SQL_LEN];

	SQLTCHAR *postgres_time_add_sql = TEXT(
		"(?::timestamp + ((? || ' SECONDS')::interval))");  
	SQLTCHAR *mysql_time_add_sql = TEXT(
		"(DATE_ADD(?, INTERVAL ? SECOND))");  

	_stprintf(postgres_time_diff_sql, TEXT(
		"extract(epoch from (time - %s))"), postgres_time_add_sql);
	_stprintf(mysql_time_diff_sql, TEXT(
		"TIME_TO_SEC(TIMEDIFF(time, %s)"), mysql_time_add_sql);  


	_stprintf(sqltext, TEXT(
		"SELECT %s, value \
		FROM %s \
		WHERE \
			cid = ? \
				AND \
			time <= %s \
		ORDER BY \
			ABS(%s) \
			ASC LIMIT 1;"), 
			postgres_time_diff_sql, 
			dat_tab,
			postgres_time_add_sql,
			postgres_time_diff_sql );

	SQLRETURN rc;
	if (!hStmt_dat_prepared) {
		rc = SQLPrepare(hStmt_dat, sqltext,  SQL_NTS);
		if (rc) {
			handleDiagnosticRecord(hStmt_dat, SQL_HANDLE_STMT, rc);
			if (rc==SQL_ERROR) return CANT_PREPARE_QUERY;
		}
		hStmt_dat_prepared = 1;
	}

	rc = SQLBindParameter(hStmt_dat, 1, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, &timestamp, 0, NULL);
	rc = SQLBindParameter(hStmt_dat, 4, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, &timestamp, 0, NULL);
	rc = SQLBindParameter(hStmt_dat, 6, SQL_PARAM_INPUT, 
		SQL_C_TIMESTAMP, SQL_TYPE_TIMESTAMP, 23, 3, &timestamp, 0, NULL);
	rc = SQLBindParameter(hStmt_dat, 2, SQL_PARAM_INPUT, 
		SQL_C_LONG, SQL_INTEGER, 10, 0, &timeshift, 0, NULL);
	rc = SQLBindParameter(hStmt_dat, 5, SQL_PARAM_INPUT, 
		SQL_C_LONG, SQL_INTEGER, 10, 0, &timeshift, 0, NULL);
	rc = SQLBindParameter(hStmt_dat, 7, SQL_PARAM_INPUT, 
		SQL_C_LONG, SQL_INTEGER, 10, 0, &timeshift, 0, NULL);

	rc = SQLBindParameter(hStmt_dat, 3, SQL_PARAM_INPUT, 
		SQL_C_LONG, SQL_INTEGER, 10, 0, &chan_key, 0, NULL);

	rc = SQLExecute(hStmt_dat);
	if (rc) {
		handleDiagnosticRecord(hStmt_dat, SQL_HANDLE_STMT, rc);
		if (rc==SQL_ERROR) return CANT_QUERY_DATA;
	}

	SQLBindCol(hStmt_dat, 1, SQL_C_LONG, time_offset, 0, NULL);
	SQLBindCol(hStmt_dat, 2, SQL_C_DOUBLE, data_out, 0, NULL);

	rc = SQLFetch(hStmt_dat);
	if (rc) {
		handleDiagnosticRecord(hStmt_dat, SQL_HANDLE_STMT, rc);
		if (rc==SQL_ERROR) return NO_DATA_FOR_THE_CHANNEL;
	}
	SQLCloseCursor(hStmt_dat);
	return OK;

}

void Provider::handleDiagnosticRecord(SQLHANDLE      hHandle,    
									  SQLSMALLINT    hType,  
									  RETCODE        RetCode)
{
	SQLSMALLINT iRec = 0;
	SQLINTEGER  iError;
	TCHAR       wszMessage[1000];
	TCHAR       wszState[SQL_SQLSTATE_SIZE+1];


	if (RetCode == SQL_INVALID_HANDLE)
	{
		_ftprintf(stderr, TEXT("Invalid handle!\n"));
		return;
	}

	while (SQLGetDiagRec(hType,
		hHandle,
		++iRec,
		wszState,
		&iError,
		wszMessage,
		(SQLSMALLINT)(sizeof(wszMessage) / sizeof(TCHAR)),
		(SQLSMALLINT *)NULL) == SQL_SUCCESS)
	{
		// Hide data truncated..
		if (wcsncmp(wszState, TEXT("01004"), 5))
		{
			_ftprintf(stderr, TEXT("[%5.5s] %s (%d)\n"), wszState, wszMessage, iError);
		}
	}


}

Provider::Err Provider::loadChannels(Network* net) {
	RETCODE rc;
	SQLTCHAR sqltext[MAX_SQL_LEN];
	_stprintf(sqltext, 
		TEXT("SELECT id, net_id_str, msmt_t, l_lim, r_lim, stde FROM %s;"), chan_tab);


	SQLUINTEGER id; 
	SQLCHAR net_id[MAX_NET_ID_LEN];
	SQLCHAR type[2];
	SQLDOUBLE llim, rlim; //left (lower) and right limits of the sensor
	SQLDOUBLE stde; // sensor std err of measurements
	SQLLEN	id_indi, netid_indi, type_indi, ll_indi, rl_indi, s_indi; //indicators

	SQLBindCol(hStmt, 1, SQL_C_ULONG, &id, 0, &id_indi);
	SQLBindCol(hStmt, 2, SQL_C_CHAR, &net_id, MAX_NET_ID_LEN, &netid_indi);
	SQLBindCol(hStmt, 3, SQL_C_CHAR, &type, 2, &type_indi);
	SQLBindCol(hStmt, 4, SQL_C_DOUBLE, &llim, 0, &ll_indi);
	SQLBindCol(hStmt, 5, SQL_C_DOUBLE, &rlim, 0, &rl_indi);
	SQLBindCol(hStmt, 6, SQL_C_DOUBLE, &stde, 0, &s_indi);

	if (rc = SQLExecDirect(hStmt, sqltext, SQL_NTS)) {
		handleDiagnosticRecord(hStmt, SQL_HANDLE_STMT, rc);
		if (rc==SQL_ERROR) return CANT_LOAD_CHANNELS;
	}

	SQLLEN n_rows = 0;
	SQLRowCount(hStmt, &n_rows);
	if (n_rows == 0) return CANT_LOAD_CHANNELS;

	nChan = (unsigned)n_rows;

	Channel* chan_list_head = NULL;

	while ((rc = SQLFetch(hStmt)) != SQL_NO_DATA) {

		Channel* aChan = new Channel();
		aChan->key = id;  
		strncpy(aChan->name, (char*)net_id, MAX_NET_ID_LEN);
		aChan->provider = this;
		switch (type[0]) {
		case 'C': aChan->type = Channel::C; break;
		case 'L': aChan->type = Channel::L; break; //actually, tank/reservior level
		case 'Q': aChan->type = Channel::Q ;break;
		case 'P': aChan->type = Channel::P ;break; 
		case 'A': aChan->type = Channel::A ;break; 
		case 'B': aChan->type = Channel::B ;break; 
		case 'D': aChan->type = Channel::D ;break; 
		case 'F': aChan->type = Channel::F ;break; 
		case 'V': aChan->type = Channel::V ;break; 
		default: aChan->type = Channel::NONE; 
		}

		aChan->lower_lim = llim;
		aChan->upper_lim = rlim;
		aChan->stde = stde;
        if (net!=NULL) {
			int idx = net->name2idx(aChan->name);
            if (idx==0) return CANT_FIND_COR_COMP;
			aChan->mindex = idx;
		} else {
            aChan->mindex = 0; //network unset
		}
		//unit type unset
		// rse, err_data should be defined 
		aChan->status = Channel::OK;

		//insert the channel
		aChan->next = chan_list_head;
		chan_list_head = aChan;
	}

	lsChan = chan_list_head;


	SQLCloseCursor(hStmt);

	return OK;

}

Provider::~Provider() {
	//destructor
	if (hStmt) { SQLFreeHandle(SQL_HANDLE_STMT, hStmt); }
	if (hStmt_dat) { SQLFreeHandle(SQL_HANDLE_STMT, hStmt); }
	if (hStmt_w) { SQLFreeHandle(SQL_HANDLE_STMT, hStmt); }

	if (hDbc)
	{
		SQLDisconnect(hDbc);
		SQLFreeHandle(SQL_HANDLE_DBC, hDbc);
	}

	if (hEnv)
	{
		SQLFreeHandle(SQL_HANDLE_ENV, hEnv);
	}


	Channel* head = lsChan; 
	Channel* cur;
	while (nChan--) {
		cur = head;
		head = head->next;
		delete cur;
	}

	//_ftprintf(stderr, TEXT("\nDisconnected."));
}


