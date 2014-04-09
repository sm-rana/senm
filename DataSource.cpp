#include <Windows.h>
#include <sqlext.h>
#include <atlconv.h>
#include "iniparser.h"
#include "DataSource.h"


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
	} else (*ds)->n_prov = nprov;

    for (int iprov=1; iprov<=nprov; ++iprov) {
        char src[MAX_DSN_LEN];
        char* dsn, *server, *port, *uid, *pwd;
        sprintf(src, "Source%d:dsn", iprov);
		dsn = iniparser_getstring(dict, src, "");

        sprintf(src, "Source%d:server", iprov);
		server = iniparser_getstring(dict, src, "");

        sprintf(src, "Source%d:port", iprov);
		port = iniparser_getstring(dict, src, "");

        sprintf(src, "Source%d:uid", iprov);
		uid = iniparser_getstring(dict, src, "");

        sprintf(src, "Source%d:pwd", iprov);
		pwd = iniparser_getstring(dict, src, "");

        sprintf(src, "DSN=%s; SERVER=%s; PORT=%s; UID=%s; PWD=%s", 
			dsn, server, port, uid, pwd);

		TCHAR tdsn[MAX_DSN_LEN];
		_tcscpy(tdsn, A2T(src));

        Provider* aPrv;
        if (Provider::Err perr = Provider::New(tdsn, &aPrv)) {
			Provider::report(perr);
            err = CANNOT_CREATE_PROVIDER;
            goto CLEAN;
		}

		if (Err derr = (*ds)->addProvider(aPrv)) {
            report(derr);
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

void DataSource::report(Err err, TCHAR* arr) {
	TCHAR errtxt[DUMMY_LAST][DATASOURCE_MAX_ERR_TEXT] = {
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

	_tcscpy(arr, errtxt[err]);
}

void DataSource::report(Err err, FILE* errfile) {
	TCHAR tperrtxt[DATASOURCE_MAX_ERR_TEXT];
	report(err, tperrtxt);
	_ftprintf(errfile, TEXT("Data Source: %s\n"), tperrtxt);
}



DataSource::Err DataSource::addProvider(Provider* prov_in) {
	unsigned n_new_chan;
    Provider::Err perr;

	if (prov_in == NULL) return PROVIDER_INVALID;
	if (perr = prov_in->loadChannels()) {
        Provider::report(perr);
		return CANT_ADD_PROVIDER;
	}
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
	TCHAR txts[10][DATASOURCE_MAX_CHAN_INFO] = {
		TEXT("No channel"),
		TEXT("Water level"),
		TEXT("Pump flow rate"),
		TEXT("Minor headloos coefficient"),
		TEXT("Discharge-side pressure"),
		TEXT("Suction-side pressure"),
		TEXT("Pipe flow rate "),
		TEXT("Pipe/valve on/off"),
		TEXT("Nodal pressure"),
		TEXT("Real-time demad (flow meter)")};

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
		default: type= TEXT("Unknown measurement");
		}

		switch (chan_it->status) {
		case Channel::OK: status = TEXT("OK"); break;
		case Channel::DATA_ERR: case Channel::DATA_OUTBOUND:
			status = TEXT("Anormaly detected"); break;
		}

		_ftprintf(stderr, TEXT(
			"Channel [%d]: Provider: %d, EPANET id: %S, Type: %s, Status: %s\n"), 
			chan_it->key, (unsigned)chan_it->provider, 
			chan_it->name, type, status);
	}
}

DataSource::Err DataSource::fillSnapshots(Tstamp t_in, unsigned n, double* snapshots) {
	//snapshots memory must be alloced,
	if (snapshots == NULL) return MEM_NOT_ALLOCED;
	unsigned j = 0;  //j - channel no.
    int fErr = 0;

	for (unsigned i=0; i<=n; ++i) { // iterate through time steps
		j=0;
		for (Channel* chan_it=lsChan; 
			chan_it; chan_it=chan_it->next) { 
				Provider::Err err = chan_it->provider->getDataAt(
					t_in, -(int)i*dt, /* time shifts backward*/
					chan_it->key, &snapshots[i*n_chan + (j++)], NULL);
				if (err) {
					Provider::report(err);
                    fErr = 1;
				}
		}
	}

    if (fErr) return CANT_GET_DATA_THRU_CHAN;
	else return OK;
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
	_tcscpy(dat_tab, TEXT("Msmts"));  // fact table
	_tcscpy(chan_tab, TEXT("Channels"));  // data table

}

Provider::Err Provider::New(TCHAR* tdsn, Provider** prov_out) {
    (*prov_out)  = new Provider();
    if (*prov_out == NULL) return MEM_NOT_ALLOCED;

    Err err;
	err = (*prov_out)->connect(tdsn);
    if (err)  goto END;

	err = (*prov_out)->loadChannels();
    if (err) goto END;

	return OK;

END:
    delete (*prov_out);
    return err;

}

void Provider::report(Err err, TCHAR* arr) {
	TCHAR errtxt[DUMMY_LAST][PROVIDER_MAX_ERRTXT] = {
		TEXT("OK"),
		TEXT("Memory can not be allocated."),

		TEXT("Unable to allocate an environment handle"),
		TEXT("can not set env handle to v3"),
		TEXT("can not allocate db connection handle"),
		TEXT("DSN string is either null or too long."),
		TEXT("Cannot connect to database. Check the DSN. "),
		TEXT("Unable to allocate statement handle"),
		TEXT("Cannot load data channels. check if table and column exist in the db."),
		TEXT("the type of the measurments is unkonwn."),
		TEXT("Database, statment, or environement handle is not ready. run connect() first."),

		TEXT("SQL statement fails to execute. Do the required columns exist in the table?"),
		TEXT("SQL statement fails to fetch"),
		TEXT("The SQL statment is not prepared. Do the data table/view exist in the db?"),
	};

	_tcscpy( arr, errtxt[err]);
}

void Provider::report(Err err, FILE* errfile) {
	TCHAR tperr[PROVIDER_MAX_ERRTXT];
	report(err, tperr);
	_ftprintf(errfile, TEXT("Provider: %s\n"), tperr);
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

	return OK;

}



Provider::Err Provider::getDataAt(SQL_TIMESTAMP_STRUCT timestamp, 
								  int timeshift,  /* desired time shift of time stamp, in sec*/
								  int chan_key,   /* the desired channel*/
								  double* data_out,    /* The data */
								  int* time_offset /*timestamp of the data returned */ 
								  ) 
{
	if ((hEnv == NULL) || (hDbc = NULL) || (hStmt == NULL) )
		return DB_CONN_NOT_READY;

	SQLTCHAR sqltext[MAX_SQL_LEN];
	_stprintf(sqltext, TEXT(
		"SELECT\
		TIME_TO_SEC(TIMEDIFF(time, DATE_ADD(?, INTERVAL ? SECOND))), \
		value \
		FROM\
		%s \
		WHERE \
		cid = ? \
		AND \
		time <= DATE_ADD(?, INTERVAL ? SECOND)\
		ORDER BY \
		ABS(TIME_TO_SEC(TIMEDIFF(\
		time, DATE_ADD(?, INTERVAL ? SECOND)))) \
		ASC LIMIT 1;"), dat_tab);

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

Provider::Err Provider::loadChannels( ) {
		RETCODE rc;
		SQLTCHAR sqltext[MAX_SQL_LEN];
		_stprintf(sqltext, 
			TEXT("SELECT id, net_id_str, msmt_t, l_lim, r_lim FROM %s;"), chan_tab);


		SQLUINTEGER id; 
		SQLCHAR net_id[MAX_NET_ID_LEN];
		SQLCHAR type[2];
		SQLDOUBLE llim, rlim; //left (lower) and right limits of the sensor
		SQLLEN	id_indi, netid_indi, type_indi, ll_indi, rl_indi; //indicators

		SQLBindCol(hStmt, 1, SQL_C_ULONG, &id, 0, &id_indi);
		SQLBindCol(hStmt, 2, SQL_C_CHAR, &net_id, MAX_NET_ID_LEN, &netid_indi);
		SQLBindCol(hStmt, 3, SQL_C_CHAR, &type, 2, &type_indi);
		SQLBindCol(hStmt, 4, SQL_C_DOUBLE, &llim, 0, &ll_indi);
		SQLBindCol(hStmt, 5, SQL_C_DOUBLE, &rlim, 0, &rl_indi);

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
			aChan->mindex = 0; //unset, deferred to network assignment
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
	if (hStmt)
	{
		SQLFreeHandle(SQL_HANDLE_STMT, hStmt);
	}

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


