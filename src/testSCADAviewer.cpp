#include <tchar.h>
#include <atltime.h>
#include <Windows.h>
#include "rvgs.h"
#include "SCADAviewer.h"


int main(int argc, char** argv) {
USES_CONVERSION;
	LPCTSTR inp = TEXT( "network.inp");
	const char* data_source_config_file = "scada_ds.ini";

	// argv[1] : inp file
	// argv[2] : ini file for database configs
	if (argc > 1) {
		inp = A2T(argv[1]);
	}

	if (argc > 2) {
		data_source_config_file = argv[2];
	}

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

	SCADAviewer *senm = new SCADAviewer();
	if (senm->init(inp, data_source_config_file, t0, 900, varima_weekly) == SCADAviewer::OK) {

		//senm->init(inp, NULL, t0, 900, varima_weekly);
		senm->run();
	}

}
