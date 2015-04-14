#include "SenmCoreIncs.h"
#include "ScadaGenerator.h"
#include "Solver.h"


// this test uses the Solver class to conduct a steady-state hydraulic simulation
// the solver is hooked into a ScadaGenerator, which synthesizes the demand and the 
// hydraulics and populates the scada database. the solver uses the generated demand 
// and get control information from the database and run hydraulic simulation. 
// The results are compared with thosed sent from the P and Q channels


void testSolver(void*, ScadaGenerator*);
Varima* makeTestVarima(int nuser, float* ref_d); 

int main(int argc, char** argv) {
	USES_CONVERSION;

	_tprintf(TEXT("Test Class Solver and Network. \n"));

	Network *testnet;  
	DataSource *ds;
	Solver* asolver;

	Network::ErrorCode ecn = Network::getNetwork(A2T(argv[1]), &testnet);
	if (ecn != Network::OK) {
		ewi(TEXT("Could not create network."));
		return 1;
	}
	testnet->report();

	DataSource::Err ecd = DataSource::New(argv[2], testnet, &ds);
	if (ecd != DataSource::OK) {
		ewi(TEXT("Could not create the data source."));
		return 2;
	}
	ds->dumpChannelsInfo();

	Solver::EWICode ecs = Solver::createSolver(testnet, ds, 0, &asolver);
	if (ecs != Solver::OK) {
		ewi(TEXT("Could not create the solver."));
		return 3;
	}

	// create a varima model
	Varima* md = makeTestVarima(testnet->getNusers(), testnet->getRefD());

	// create a scada generator and set hook function
	ScadaGenerator* sg;
	SG_ERR ecsg = SG_new(A2T(argv[1]), argv[2], 1234, md, &sg);
	sg->extObj = (void*)asolver;
	sg->extUpdate = testSolver;

	//Run generator
	Tstamp t0 = { 2013, 4, 22,  0, 0, 0, 0 };
	int timespan = 3600*24*7;
	SG_make(sg, &t0, timespan); 

	SG_delete(sg);

	return 0;
}

void testSolver(void* sol_in, ScadaGenerator* sg) {
	Solver* solver = (Solver*)sol_in;

	// update data source buffer
	solver->_ds->updateBufSnapshot(sg->t0 + CTimeSpan(0,0,0,sg->elapTime));

	double* XD = (double*)calloc(solver->_nXd, sizeof(double));
	int ijunc; int iuser; double tpD;
	for (ijunc=1, iuser=0; ijunc<=sg->njunc; ++ijunc) {
		if (sg->demands[ijunc])  tpD = sg->_tpd[iuser++]; // get demand data 
		int ixd = solver->_tabIdx2Xd[ijunc]; // fill the XD array (so no D channel)
		if (ixd != -1) XD[ixd] = tpD;
	}

	// TODO: Fix the solver test
	solver->run(XD, solver->_nXd);
    free(XD);

	// compare modeling result and P, D channels
	int iichan = 0;
	for (Channel* ichan=solver->_ds->lsChan; ichan!=NULL; ichan=ichan->next, ++iichan) {
		if (ichan->type == Channel::P)  {
			printf("P at Node %d: EPANET %6.2f ,  SenM Solver %6.2f \n", ichan->mindex,
				solver->_ss[iichan]/solver->_net->Ucf[Network::PRESSURE]+ 
				solver->_net->Node[ichan->mindex].El, solver->H[ichan->mindex]);
		}
		if (ichan->type == Channel::Q) {
			printf("Q at Node %d: EPANET %6.2f, SenM Solver %6.2f \n", ichan->mindex,
				solver->_ss[iichan]/solver->_net->Ucf[Network::FLOW], 
				solver->Q[ichan->mindex]);
		}
	}

	//Sleep(5000);
}




Varima* makeTestVarima(int nuser, float* ref_d) {
	// the model structure described here comes from research of (chen 2014)
	//the demand model has 15-min interval 
	// input/output of the model is log(demand)

	Varima* vsar_h = new Varima(nuser, 1, 0, 0);  //p, q, d
	Varima* vsar_d = new Varima(vsar_h, 1, 0, 24*4, 1); //P,Q,S,D
	Varima* vsar_w = new Varima(vsar_d, 1, 0, 168*4, 1);

	double cv = 0.0743;  //coefficient of  variation
	double* phi_h = new double[nuser * 1];
	double* phi_d = new double[nuser * 1];
	double* phi_w = new double[nuser * 1];

	double* cov = new double[nuser*nuser];

	//random number generator for covariance, don't know appropriate values
	Rvgs* rgp = new Rvgs(1111); 

	for (int i=0; i<nuser; ++i) {
		phi_h[i] = 0.6001;
		phi_d[i] = -0.4582;

		phi_w[i] = -0.7651; //  \Phi_1

		for (int j=i; j<nuser; ++j) {
			if (j==i) cov[i*nuser + j] = cov[j*nuser+i] = SQR(cv * ref_d[i]);
			else { // set covariances randomly 
				cov[i*nuser + j]=cov[j*nuser+i] = 
					cv * ref_d[i] * ref_d[j] * rgp->Uniform(-0.5, 0.5);
			}
		}
	}

	// test varima and generator
	vsar_h->setPara(phi_h, NULL, cov);
	vsar_d->setPara(phi_d, NULL, NULL);
	vsar_w->setPara(phi_w, NULL, NULL);

	//return vsar_w;
	return vsar_d;

}