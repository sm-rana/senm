/*Test ScadaGenerator.cpp, 
run sql_clean.sql and ctown_channel_c.sql before running this test*/


#include "ScadaGenerator.h"
#include "f2c.h"
#include "clapack.h"

int main () {

	LPCTSTR inp = TEXT("Data\\c-town_true_network.inp");
	LPCTSTR dsn = TEXT(
		 "DSN=jcx;\
		 DESCRIPTION={Jinduan's Senm Database};\
		 SERVER=localhost;\
		 UID=rtx_db_agent;\
		 PWD=rtx_db_agent;\
		 DATABASE=rtx_demo_db;\
		 PORT=3306;\
		 FOUND_ROWS=1");


	//// make iid lognormal demands
	//ScadaGenerator* psg = new ScadaGenerator();
	//psg->init(inp, dsn, 1234);

	//SQL_TIMESTAMP_STRUCT t0 = {
	//	2013, 4, 22,  
	//	4, 7, 0, 0
	//};

	//psg->make(&t0, 3600*24*14 , 0.03);
	//delete psg;
	

	// make spatial-temporal demands
	ScadaGenerator* psg2 = new ScadaGenerator();
	psg2->init(inp, dsn, 1234);

	int nuser = psg2->getNuser(); //388 in this example

	Varima* varima = new Varima(nuser, 2, 0);  //p, q, d
	//the network has 15-min demand interval so daily period = 24*4
	Varima* varima_diurnal = new Varima(varima, 1, 0, 24*4, 1); //P,Q,S,D
	Varima* varima_weekly = new Varima(varima_diurnal, 1, 0, 24*4*7, 1);

	double sigma = 0.0003;  // variance of log(demand)
	double* phi_in = new double[nuser*2];
	double* phi_in2 = new double[nuser];
	double* phi_in3 = new double[nuser];
	double* cov = new double[nuser*nuser];

	Rvgs* rgp = new Rvgs(1111); //random number generator for parameters
	for (int i=0; i<nuser; ++i) {
		switch (i%3) {
		case 0: 
			phi_in[2*i] = 1.28;
			phi_in[2*i+1] = -0.4;
			break;
		case 1:
			phi_in[2*i] = 0.7;
			phi_in[2*i+1] = 0;
			break;
		case 2:
			phi_in[2*i] = 1;
			phi_in[2*i+1] = -0.5;
			break;
		}

		switch (i%2) {
		case 0:
			phi_in2[i]=-0.2;
			break;
		case 1:
			phi_in2[i]=0.12;
			break;
		}

		switch (i%5) {
		case 0: phi_in3[i] = 0.2; break;
		case 1: phi_in3[i] = -0.3; break;
		case 2: phi_in3[i] = rgp->Normal(0.1, 0.1); break;
		case 3: phi_in3[i] = -0.14; break;
		case 4: phi_in3[i] = rgp->Normal(0, 0.2); break;
		}

		for (int j=i; j<nuser; ++j) {
			if (j==i) cov[i*nuser + j] = cov[j*nuser+i] = sigma;
			else cov[i*nuser + j]=cov[j*nuser+i]=sigma*rgp->Uniform(-0.5, 0.5);
		}
	}

	//test gotoblas
	//double alpha = 1; integer incx=2; double beta=0; integer incy=1;
	//double* y = new double[nuser];
	//dgemv_("N", (integer*)&nuser, (integer*)&nuser, 
	//	   &alpha, 
	//	   cov, (integer*)&nuser, 
	//	   phi_in, &incx, 
	//	   &beta, 
	//	   y, &incy); // U matrix * vector

	// test varima and generator
	varima->setPara(phi_in, NULL, cov);
	varima_diurnal->setPara(phi_in2, NULL, NULL);
	varima_weekly->setPara(phi_in3, NULL, NULL);
	
	
	SQL_TIMESTAMP_STRUCT t0_2 = {
		2013, 	5, 		29,
		4,  7, 	0, 	0
	};
	psg2->make(&t0_2, 3600*24*14 /*2 weeks*/, varima_weekly/*cov*/);

	delete psg2; delete rgp;
	delete varima; delete varima_diurnal; delete varima_weekly;


}
