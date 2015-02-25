// EMMC program of senm system

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "iniparser.h"
#include "rvgs.h"

#include "SenMCoreIncs.h"

#include "Network.h"
#include "DataSource.h"
#include "Solver.h"
#include "Population.h"
#include "VarModel.h"


int main(int argc, char** argv) {
	USES_CONVERSION;

	if (argc != 3) return 1;

	int mode = 0; // data integrity check
	if (strcmp(argv[1], "-estimate") == 0) {
		mode = 1; // emmc
	}

	printf("[0]Loading config file %s ... \n", argv[2]);
	dictionary* dict;
	dict = iniparser_load(argv[2]);
	if (dict == NULL) {
		return 2;
	}

	char* inp_file_ascii = iniparser_getstring(dict, "Net:inp_file", "");
	if (strcmp(inp_file_ascii, "") ==0) {
		return 3;
	}
	TCHAR inp_file_tchar[MAX_FILE_NAME_SIZE];
	_tcscpy(inp_file_tchar, A2T(inp_file_ascii));

	char* db_config_file = iniparser_getstring(dict, "Scada:db_conn_file", "");
	if (strcmp(db_config_file, "") ==0) {
		return 5;
	}

	char* xd_model_file = iniparser_getstring(dict, "XDemands:model_file", "");
	if (strcmp(xd_model_file, "") ==0) {
		return 5;
	}

	int mcmc_chain_size = iniparser_getint(dict, "EMMC:max_chain_size", 0);
	if (mcmc_chain_size <= 0) {
		return 101;
	}

	int est_win_size = iniparser_getint(dict, 
		"EMMC:estimation_data_window_size", 0);
	if (est_win_size <= 0) {
		return 102;
	}

	double prop_std = iniparser_getdouble(dict, 
		"EMMC:proposal_std", 0);
	if (prop_std <=0) {
		return 103;
	}

	char* est_end_time = iniparser_getstring(dict, 
		"EMMC:estimation_end_time", "");
	if (strcmp(est_end_time, "") == 0) {
		return 13;
	}

	
	Network* net;
	DataSource* ds;
	Solver* sv;
	Population* pop;
	VarModel* vm;


	printf("[0]Loading network file %s ... \n", inp_file_ascii);
	Network::ErrorCode nec = Network::getNetwork(inp_file_tchar, &net);
	if (nec != Network::OK) {
		return 4;
	}

	if (mode == 0) net->report();

	printf("[0]Establishing Scada DB connection with config file %s ...\n",
			db_config_file);
	DataSource::Err dsec = DataSource::New(db_config_file, net, &ds);
	if (dsec != DataSource::OK) {
		return 6;
	}


	printf("[0]Building hydraulic simulation solver ...\n");
	Solver::EWICode sec = Solver::createSolver(net, ds, 0, &sv);
	if (sec != Solver::OK) {
		return 7;
	}

	printf("[0]Building a xdemand model specified in file %s...\n",
			xd_model_file);
	int dim;
	int ns;
	int* s;
	int* p;
	int np;
	double *mu;
	double *phi0, *phi; // row- and column- major
	double *cov0, *cov;
	FILE* xmf;

	xmf = fopen(xd_model_file, "r");
	fscanf(xmf, "%d", &dim);
	fscanf(xmf, "%d", &ns);

	if (dim<=0 || ns <=0) return 8;

	s = (int*)malloc(sizeof(int)*(ns+1));
	p = (int*)malloc(sizeof(int)*(ns+1));

	for (int is=0; is<=ns; ++is) {
		fscanf(xmf, "%d", &s[is]);
	}

	np = 0;
	for (int is=0; is<=ns; ++is) {
		fscanf(xmf, "%d", &p[is]);
		np += p[is];
	}

	mu = (double*)calloc(dim, sizeof(double));
	cov0 = (double*)calloc(dim*dim, sizeof(double));
	cov = (double*)calloc(dim*dim, sizeof(double));
	phi0 = (double*)calloc(dim*dim*np, sizeof(double));
	phi = (double*)calloc(dim*dim*np, sizeof(double));

	for (int i = 0; i < dim; ++i ) {
		fscanf(xmf, "%lf", &mu[i]);
	}
	for (int i=0; i<dim*dim*np; ++i) {
		fscanf(xmf, " %lf", &phi0[i]);
	}

	for (int i =0; i<dim*dim; ++i) {
		fscanf(xmf, " %lf", &cov0[i]);
	}

	double *xd0; //initial water demand panel
	xd0 = (double*)calloc(dim*est_win_size, sizeof(double));
	for (int i=0; i<dim*est_win_size; ++i) {
		fscanf(xmf, "%lf", &xd0[i]);
	}

	fclose(xmf);

	// row major to column major
	for (int ip=0; ip<np; ++ip) 
		for (int ir =0; ir<dim; ++ir) 
			for (int ic = 0; ic<dim; ++ic) {
				phi[ip*dim*dim + ic*dim + ir] = 
					phi0[ip*dim*dim + ir*dim + ic];
			}

	for (int ir =0; ir<dim; ++ir) 
		for (int ic = 0; ic<dim; ++ic) {
			cov[ic*dim + ir] = cov0[ir*dim + ic];
		}

	VarModel_Err vmec = VarModel_new(dim, ns, s, p, &vm);
	VarModel_setMu(vm, mu, dim, 1);
	VarModel_setPhi(vm, phi, dim*dim*np);
	VarModel_Err vmec_cov = VarModel_setCov(vm, cov, dim*dim);

	if (mode == 0) VarModel_dump(vm);

	free(s); free(p);
	free(phi); free(phi0); free(cov0); free(cov);

	if (vm->n != sv->_nXd) {
		return 10;
	}

	if (mode == 0) {
		printf("  EPANET id (idx) <-> Element idx of VAR model mapping\n");
		for (int ixd = 0; ixd < sv->_nXd; ++ixd) {
			char eid[MAX_NET_ID_LEN];
			int eidx = sv->_tabXd[ixd];
			net->nodeId2NodeName(eidx, eid, MAX_NET_ID_LEN);
			printf(" %s (%d) <-> %d\n", eid, eidx, ixd);
		}
		printf("\n");
	}


	printf("[0]Building water demand population ...\n");
	Pop_Err pec = Pop_new(mcmc_chain_size, dim, est_win_size, &pop);

	if (mode == 0) {
		printf("SenM checked out successfully.");
		goto END;
	}

	// Loading scada panel data

	printf("[0]Loading scada panel data from the channels ...\n");
	printf("  Time interval %d s, %d snapshots till (%s) used in EM. \n", 
		ds->dt, est_win_size, est_end_time);
	Tstamp tcur;  // estimation end time for EM, usually the current time stamp
	sscanf(est_end_time, "%d %d %d %d %d %d", 
		&tcur.year, &tcur.month, &tcur.day, 
		&tcur.hour, &tcur.minute, &tcur.second);
	
	double* scada_panel;
	int n_ch = ds->n_chan;
	scada_panel = (double*)calloc(n_ch*est_win_size, sizeof(double));
	dsec = ds->fillSnapshots(tcur, est_win_size, scada_panel);


	//first sample
	int ws = est_win_size;
	memcpy(&POPH2(pop, 0, 0, 0, 0), xd0, sizeof(double)*dim*ws);
	free(xd0);

	//init random number generator
	Rvgs* rng = new Rvgs(1);

	// MCMC
	double cur_ll = -DBL_MAX;  //current log-likelihood
	double last_ll = -DBL_MAX;  //ll for the previous sample
	for (int im = 0; im < mcmc_chain_size - 1; ++im) {

		cur_ll = 0;
		//compute hydraulic model likelihood
		for (int it=0; it<ws; ++it) { //each time step in the est window
			double* cur_xd = &POPH2(pop, 0, im, 0, it);
			double* cur_ch_data = &scada_panel[it*n_ch];
			sv->run(cur_xd, dim, cur_ch_data, n_ch);
			cur_ll += sv->logL(cur_ch_data, n_ch);
		}

		//compute demand model likelihood
		double dm_ll=0;
		int nea = 0;
		double* cur_xd_panel = &POPH(pop, 0, im, 0);
		VarModel_logLikelihood(vm, cur_xd_panel, dim, ws, &dm_ll, &nea);
		cur_ll += dm_ll;

		//MH accept-reject
		if (cur_ll < last_ll && rng->Random() > exp(cur_ll - last_ll)) {
			// reject, roll-back to previous sample
			memcpy(cur_xd_panel, &POPH(pop, 0, im-1, 0), sizeof(double)*dim*ws);
		}

		// report chain stats
		if ((im+1) % 200 == 0) Pop_report(pop);

		last_ll = cur_ll;
		//propose a new demand panel
		for (int id=0; id<dim*ws; ++id) {
			POPH(pop, 0, im+1, id) = POPH(pop, 0, im, id) 
				+ rng->Normal(0, prop_std);
		}

	}

	Pop_report(pop);



	free(scada_panel);
	delete rng;

END:
	// clean up

	VarModel_del(&vm);
	Pop_del(&pop);
	delete sv;
	delete ds;
	delete net;

	iniparser_freedict(dict);

	return 0;


}


