// Test Scada Generator 
// If vis = 1, Use Vtk to visualize network demands and hydraulics
#include <process.h>
#include <Windows.h>
#include "Varima.h"
#include "ScadaGenerator.h"
#include "VisScadaGenerator.h"

int vis = 0;

// get example varima model
Varima* makeTestVarima(int, float*);

int main(int argc, char** argv) {
USES_CONVERSION;
    //argv[1]: inp file
    //argv[2]: ini file
	//argv[3]: "visual" generator

	if (argc >= 4) {
		if (strcmp(argv[3], "-v") == 0 || strcmp(argv[3], "-V") == 0 ) {
			vis = 1;
		} else {
			vis = 0;
		}
	}


    //start time and time span
	Tstamp t0 = { 2013, 4, 22,  0, 0, 0, 0 };
	int timespan = 3600*24*30;

    //read inp file (1st time for nusers)
	Network::ErrorCode ec;
    Network *net;
	ec = Network::getNetwork(A2T(argv[1]), &net);
	Network::report();
    if (ec != Network::OK) return 1;

	// create a varima model
	Varima* varima_weekly = makeTestVarima(net->getNusers(), net->getRefD());

	if (vis == 1) {
		// create Visualizer of scada generator
		VSG_ERR vsgerr;
		VisScadaGenerator* vsg;
		vsgerr = VSG_new(A2T(argv[1]), argv[2], 1234, varima_weekly, timespan, &vsg);
		//vsgerr = VSG_new(A2T(argv[1]), argv[2], 1234, , &vsg);
		if (vsgerr!=VSG_OK) { VSG_ewi(vsgerr); return 2;}

		// run simulation and visualization
		VSG_run(vsg, &t0);

		// clean up
		VSG_delete(vsg);
	} else { // no visualization
		SG_ERR sg_err;
		ScadaGenerator* sg;
		sg_err = SG_new(A2T(argv[1]), argv[2], 1234, varima_weekly, &sg);  
		if (sg_err != SG_OK) {
			SG_ewi(sg_err); return 3;
		}
		SG_make(sg, &t0, timespan);

		SG_delete(sg);
	}

	//TODO: network destructor
    //delete net;
	return 0;

}

Varima* makeTestVarima(int nuser, float* ref_d) {
    // the model structure described here comes from research of (chen 2014)
	//the demand model has 15-min interval 
    // input/output of the model is log(demand)

	Varima* vsar_h = new Varima(nuser, 1, 0, 0);  //p, q, d
	Varima* vsar_d = new Varima(vsar_h, 1, 0, 24, 1); //P,Q,S,D
	//Varima* vsar_d = new Varima(vsar_h, 1, 0, 24*4, 1); //P,Q,S,D
	//Varima* vsar_w = new Varima(vsar_d, 1, 0, 168*4, 1);

	double cv = 0.0743;  //coefficient of  variation
	double* phi_h = new double[nuser * 1];
	double* phi_d = new double[nuser * 1];
	//double* phi_w = new double[nuser * 1];

	double* cov = new double[nuser*nuser];

    //random number generator for covariance, don't know appropriate values
	Rvgs* rgp = new Rvgs(1111); 

	for (int i=0; i<nuser; ++i) {
        phi_h[i] = 0.6001;
        //phi_h[i+nuser] = 0.76001;
        phi_d[i] = -0.4582;

		//phi_w[i] = -0.7651; //  \Phi_1
	      
		for (int j=i; j<nuser; ++j) {
			if (j==i) {
				cov[i*nuser + j] = cov[j*nuser+i] = SQR(cv * ref_d[i]);
			}
			else { // set covariances (randomly)
				cov[i*nuser + j]=cov[j*nuser+i] = 
				  cv *cv * ref_d[i] * ref_d[j] * rgp->Uniform(-0.1, 0.1);
			}
		}
	}

	// test varima and generator
	vsar_h->setPara(phi_h, NULL, cov);
	vsar_d->setPara(phi_d, NULL, NULL);
	//vsar_w->setPara(phi_w, NULL, NULL);

    //return vsar_w;
    return vsar_d;

}
