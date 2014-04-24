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

/* TestDataSource  ->  Test the class of DataSource . It is assumed that
the scada db schema (Channels and Msmts tables ) has been established in 
the database. otherwise run sql_clean.sql and c_town_channel_c.sql first.
Also the mysql service shoud be running at the time of test.*/

#include "DataSource.h"

int main(int argv, char** argc) { //argc[1] is the config file name 

    DataSource::Err derr;

    DataSource *ds;
    derr = DataSource::New(argc[1], &ds);
    if (derr) {
        DataSource::reportEWI(derr);
        return 1;
	}

	ds->dumpChannelsInfo();
    
    CTime ct = CTime(2013, 4, 22, 10, 11, 38);
	ds->updateBufSnapshot(ct);

    for (int iCh = 1; iCh<=ds->n_chan; ++ iCh) {
		printf("%d: %6.2f   ", iCh, ds->_datBuf[iCh-1]);
	}

    delete ds;
	//printf("\nTested.\n");
    return 0;


}
      

   


