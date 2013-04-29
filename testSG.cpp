/*Test ScadaGenerator.cpp, 
run sql_clean.sql and ctown_channel_c.sql first */
#include "ScadaGenerator.h"


int main () {

	const char* inp = "Data\\c-town_true_network.inp";
	const char* dsn = 
		 "DSN=jcx;\
		 DESCRIPTION={Jinduan's Senm Database};\
		 SERVER=localhost;\
		 UID=rtx_db_agent;\
		 PWD=rtx_db_agent;\
		 DATABASE=rtx_demo_db;\
		 PORT=3306;\
		 FOUND_ROWS=1";

	ScadaGenerator* psg = new ScadaGenerator();

	psg->init(inp, dsn);

	SQL_TIMESTAMP_STRUCT t0 = {
		2013,
		4,
		22,
		4,
		7,
		0,
		0
	};


	psg->make(&t0, 3600*24*14, 0.03);

	delete psg;

}
