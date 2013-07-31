#include <tchar.h>
#include <atltime.h>
#include "SenM.h"


int main() {
	LPCTSTR inp = TEXT("Data\\c-town_true_network.inp");
	LPCTSTR dsn = TEXT(
	"DSN=jcx; SERVER=localhost; UID=rtx_db_agent;PWD=rtx_db_agent; DATABASE=rtx_demo_db; PORT=3306; FOUND_ROWS=1");

	// create a test seasonal varima model
	int nuser = 334;
	double variance = 0.03;  // variance of log(demand)
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
		/*case 0: case 1: */case 2:
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
			if (j==i) cov[i*nuser + j] = cov[j*nuser+i] = variance;
			else 
				cov[i*nuser + j]=cov[j*nuser+i] = 
				  variance * rgp->Uniform(-1, 1);
		}
	}

	Varima* varima = new Varima(nuser, 2, 0);  //p, q, d
	//the network has 15-min demand interval so daily period = 24*4
	Varima* varima_diurnal = new Varima(varima, 1, 0, 24*4, 1); //P,Q,S,D
	Varima* varima_weekly = new Varima(varima_diurnal, 1, 0, 24*4*7, 1);

	// test varima and generator
	varima->setPara(phi_in, NULL, cov);
	varima_diurnal->setPara(phi_in2, NULL, NULL);
	varima_weekly->setPara(phi_in3, NULL, NULL);


	CTime* t0 = new CTime( //senm start time
	2013, 4, 22,  //year month day
	4, 7, 0); //hour minutes sec

	SenM *senm = new SenM();
	if (senm->init(inp, dsn, t0, 900, varima_weekly) == SenM::OK) {

		//senm->init(inp, NULL, t0, 900, varima_weekly);
		senm->run();
	}

}
