#include "DataSource.h"
#include <Windows.h>
#include <strsafe.h>
#include <sqlext.h>
#include "buildControl.h"

#ifdef TEST_DATA_SOURCE

int wmain() {
	LPCTSTR connection_string = TEXT("DSN=jcx;\
		 DESCRIPTION={Jinduan's Senm Database};\
		 SERVER=localhost;\
		 UID=rtx_db_agent;\
		 PWD=rtx_db_agent;\
		 DATABASE=rtx_demo_db;\
		 PORT=3306;\
		 FOUND_ROWS=1");

	Provider* provider = new Provider;
	provider->connect(connection_string);

	DataSource* ds = new DataSource(900, provider);

	ds->dumpChannelsInfo();

	unsigned sb = 3;

	double* ss = new double[ds->getNChan() * (sb + 1)];

	

	SQL_TIMESTAMP_STRUCT ts = {
		2013,
		4,
		22,
		10,   //h
		11,  //m
		38,	 //s
		0
	};


	ds->getSnapshots(ts, sb, ss);

	delete ds;
	return 0;


}

#endif