#include <Windows.h>
#include "Solver.h"
#include "Population.h"

const double Solver::MISSING = -1.0E-10;
const double Solver::QZERO = 1.e-6;

Solver::ErrorCode Solver::createSolver(Network* net, Solver** outsolver) {
	// allocate 
	
	if (net == NULL) {
		outsolver = NULL;
		return(NETWORK_NOT_EXIST);
	}

	Solver *prec = new Solver(); //solver precursor
	prec->_net = net;

	prec->N = net->MaxNodes;
	prec->M = net->MaxLinks;

	int n = prec->N + 1;
	int m = prec->M + 1;

	prec->D  = (double *) calloc(n, sizeof(double));
	prec->H =  (double *) calloc(n, sizeof(double));

	prec->Aii = (double *) calloc(n,sizeof(double));
	prec->Aij = (double *) calloc(net->Ncoeffs+1,sizeof(double));
	prec->F   = (double *) calloc(n,sizeof(double));
	prec->E   = (double *) calloc(n,sizeof(double));
	prec->P   = (double *) calloc(m,sizeof(double));
	prec->Y   = (double *) calloc(m,sizeof(double));
	prec->X   = (double *) calloc(m>n?m:n,sizeof(double));


	prec->Q    = (double *) calloc(m, sizeof(double));
    prec->K    = (double *) calloc(m, sizeof(double));
    prec->S    = (char  *) calloc(m, sizeof(char));
	

	if (prec->D==NULL || prec->H==NULL || prec->Q==NULL || 
		prec->K==NULL || prec->S==NULL) {
		delete(prec);
			outsolver = NULL;
		return(MALLOC_ERROR);
	}

/* Set exponent in head loss equation and adjust flow-resistance tolerance.*/
   if (net->Formflag == Network::HW) prec->Hexp = 1.852;
   else                prec->Hexp = 2.0;

   *outsolver = prec;

}

Solver::Solver()  //setdefaults() in input1.c
{
   Htol = 0.005;   /* Default head tolerance         */
   Hacc      = 0.001;            /* Default hydraulic accuracy     */
   MaxIter   = 200;         /* Default max. hydraulic trials  */
   ExtraIter = -1;              /* Stop if network unbalanced     */
   RQtol     = 1E-7;           /* Default hydraulics parameters  */
   CheckFreq = 2;
   MaxCheck  = 10;
   DampLimit = 0; 
}


Solver::~Solver()
{
}

Solver::ErrorCode Solver::setSnapshot(double* shot) {
	// Use the info in the (data) snapshot to set 
	// (1) link open/close status
	// (2) tank/reservoir level
	// (3) water demands at known points
	return OK;
}

Solver::ErrorCode Solver::solveHyd(double* demand) {
	int i = 0;
	ErrorCode ec = OK;

	if (demand != NULL) { // update demands



	} else { 
	// use default values of demands in the network
		for (i=1; i<=N; ++i) 
			if (D[i] == 0)  // not updated by snapshot
				if (_net->Node[i].D)  // has demand in the network junction
					D[i] = _net->Node[i].D->Base;
	}

	//initialize flows
	for (i=1; i<=M; i++)  {   /* Initialize flows */
		if (_net->Link[i].Stat == Network::CLOSED) Q[i] = QZERO;
		else if (_net->Link[i].Type == Network::PUMP) 
			Q[i] = _net->Link[i].Kc * 
			_net->Pump[(int)_net->Link[i].Diam].Q0; // set pump initial flow
		else Q[i] = 
			Network::_PI * (_net->Link[i].Diam * _net->Link[i].Diam)/4.0;  
	}

	//set tank levels


	return ec;
}

#ifdef TEST_SOLVER

int main() {

	std::cout << "Test Class Solver and Network. \n";
	Network *testnet;  Solver* asolver;
	Network::ErrorCode ec = Network::OK;
	Solver::ErrorCode ecs = Solver::OK;

	if (ec = Network::getNetwork("test\\Net3.inp", testnet)) {
		std::cout << "Can not create the network: " << ec << "\n";
		return 1;
	
	if (ecs = Solver::createSolver(testnet, asolver)) {
		std::cout << "Cannot create the solver: " << ecs << "\n";
		return 2;
	}

	ecs = asolver->setSnapshot(NULL);
	if (ecs != Solver::OK)  
		std::cout << "set controls error: " << ecs << "\n";

	ecs = asolver->solveHyd(NULL);
	if (ecs != Solver::OK)  
		std::cout << "solving hydraulics error: " << ecs << "\n";





}
#endif