// EMMC program of senm system
// jinduan:
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include <iostream>
#include <time.h>

#include "iniparser.h"
#include "rvgs.h"

#include "SenMCoreIncs.h"

#include "Network.h"
#include "DataSource.h"
#include "Solver.h"
#include "Population.h"
#include "VarModel.h"

//#include "clapack_3.2.1\clapack.h"

// Started working on it again to run the code for the modified Net3 example for 
// the journal article we are going to coauthor with Jinduan
// 3-15-2017 (6:41pm) -- SM Masud Rana

/*
 * Log of list of new edits (starting from 6-5-2017)
 * 6-5-2017 (7pm) : Obtaining random number generator seed from the ini file as a user input. Previously a fixed value of 2 was used.
 */

using namespace std;

double pdfNorm(double xx, double mean_xx, double sigma_xx)
{
	double pdf = -1;
	pdf = 1 / (sigma_xx*pow((2 * _PI), 0.5)) * exp(-0.5*pow(((xx - mean_xx) / sigma_xx), 2));
	return pdf;
}


int main(int argc, char** argv) {
	// argv[1] should be "-estimate" for mode 1
	// argv[2] should be the config file i.e., senm.ini

	USES_CONVERSION;
	
	//float senm_version = 2.00.1;	// Masud: this is just for printing out the version number in a report file. 
	
	if (argc != 3) { cout << "Missing Input Parameters. Rerun the program with two parameters.\n" << endl; return 1; }
	printf("argc = %d \n argv = %s, %s\n--\n", argc, argv[1], argv[2]);

	int mode = 0; // data integrity check
	if (strcmp(argv[1], "-estimate") == 0) {	//smmr: strcmp() returns 0 when strings exactly match!
		mode = 1; // emmc
	}
	
	if (strcmp(argv[1], "-Estep-only") == 0) {
		mode = 2; // will skip the time-series model estimation part. (added 6-27-2017)
	}
	

	printf("[0]Loading config file %s ... \n", argv[2]);
	dictionary* dict;
	dict = iniparser_load(argv[2]);
	if (dict == NULL) {
		return 2;
	}

	char* inp_file_ascii = iniparser_getstring(dict, "Net:inp_file", "");
	if (strcmp(inp_file_ascii, "") == 0) {
		return 3;
	}

	TCHAR inp_file_tchar[MAX_FILE_NAME_SIZE];
	_tcscpy(inp_file_tchar, A2T(inp_file_ascii));

	char* db_config_file = iniparser_getstring(dict, "Scada:db_conn_file", "");
	if (strcmp(db_config_file, "") == 0) {
		return 5;
	}

	char* xd_model_file = iniparser_getstring(dict, "XDemands:model_file", "");
	if (strcmp(xd_model_file, "") == 0) {
		return 5;
	}

	int is_log_demand = iniparser_getint(dict, "XDemands:is_log_demand", 0);
	int is_dev_demand = iniparser_getint(dict, "XDemands:is_dev_demand", 0);		//Masud: Not used anymore. Need to throw away.
	int is_cluster = iniparser_getint(dict, "XDemands:is_cluster", 0);

	char* actual_demand_file = iniparser_getstring(dict, "XDemands:lookback_demands", "");
	if (strcmp(actual_demand_file, "") == 0) {
		printf("\n\nNo \"actual demand\" file name specified in the MOD file.\n");
		return 1002;
	}

	char* cluster_list_file;
	if (is_cluster)
	{
		cluster_list_file = iniparser_getstring(dict, "XDemands:cluster_list", "");
		if (strcmp(cluster_list_file, "") == 0) {
			printf("\n\nNo \"cluster list\" file name specified in the MOD file.\n");
			return 1003;
		}
	}
	
	int rng_seed = iniparser_getint(dict, "EMMC:rng_seed", -1);
	if (rng_seed <= 0) {
		printf("\n\n Error reading random number generator seed from the ini file.\n\n");
		return 106;
	}
	
	int n_iter = iniparser_getint(dict, "EMMC:em_iterations", 0);
	if (n_iter <= 0) {
		return 106;
	}

	int mcmc_chain_size = iniparser_getint(dict, "EMMC:max_chain_size", 0);
	if (mcmc_chain_size <= 0) {
		return 101;
	}

	int burn_in = iniparser_getint(dict, "EMMC:burn_in", 0);
	if (burn_in < 0 || burn_in > mcmc_chain_size) {
		return 107;
	}

	int est_win_size = iniparser_getint(dict,
		"EMMC:estimation_win_size", 0); // est_win_size = 168 for Net1
	if (est_win_size <= 0) {
		return 104;
	}

	double prop_std = iniparser_getdouble(dict,
		"EMMC:proposal_std", 0);
	if (prop_std <= 0) {
		return 103;
	}

	char* est_end_time = iniparser_getstring(dict,
		"EMMC:estimation_end_time", "");
	if (strcmp(est_end_time, "") == 0) {
		return 13;
	}

	char* estep_outfile = iniparser_getstring(dict, "EMMC:estep_outfile", "");
	if (strcmp(estep_outfile, "") == 0) return 14;

	char* os_outfile = iniparser_getstring(dict, "EMMC:os_debug_outfile", "");
	if (strcmp(os_outfile, "") == 0) os_outfile = NULL;


	double norm_std_ts = iniparser_getdouble(dict, "EMMC:normdist_std", -0.5);
	if (norm_std_ts < 0)
	{
		printf("EMMC: normdist_std is not set in the MOD file.(%f)\n ", norm_std_ts);
		return 1001;
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

	if (mode == 0) net->report();	// mode = 0 for data integrity check

	printf("[0]Establishing Scada DB connection with config file %s ...\n",
		db_config_file);
	DataSource::Err dsec = DataSource::New(db_config_file, net, &ds);
	if (dsec != DataSource::OK) {
		return 6;
	}
	ds->dumpChannelsInfo();

	printf("[0]Building hydraulic simulation solver ...\n");
	Solver::EWICode sec = Solver::createSolver(net, ds, 0, &sv);
	if (sec != Solver::OK) {
		return 7;
	}

	printf("[0]Building a xdemand model specified in file %s...\n",
		xd_model_file);
	int dim; // =7 sm. number of demand variables
	int ns;	// =1 sm. number of seasons
	int* s; // = s=[1,24]
	int* p; // p=[2,1]
	int np;	// =3 sm. total number of AR params
	double *mu;
	double *phi0, *phi; // row- and column- major
	double *cov0, *cov;
	FILE* xmf;

	xmf = fopen(xd_model_file, "r");
	fscanf(xmf, "%d", &dim);	//sm./ dim == first line of the .mod file. = number of demand var
	fscanf(xmf, "%d", &ns);		//sm./ dim == second line of the .mod file. = number of seasons

	if (dim <= 0 || ns <= 0) return 8;

	s = (int*)malloc(sizeof(int)*(ns + 1));	// ns = 1 for Net1
	p = (int*)malloc(sizeof(int)*(ns + 1));

	//sm./ ns = 1 for NET1
	for (int is = 0; is <= ns; ++is) {
		fscanf(xmf, "%d", &s[is]);
	}

	np = 0;
	for (int is = 0; is <= ns; ++is) {
		fscanf(xmf, "%d", &p[is]);
		np += p[is];
	}//sm./ total number of AR parameters

	mu = (double*)calloc(dim, sizeof(double)); //sm./ dim = 7 for Net1
	cov0 = (double*)calloc(dim*dim, sizeof(double));
	cov = (double*)calloc(dim*dim, sizeof(double));
	phi0 = (double*)calloc(dim*dim*np, sizeof(double));
	phi = (double*)calloc(dim*dim*np, sizeof(double));

	FILE *fp_Spit; // "spitOut.txt" 
	fopen_s(&fp_Spit, "spitOut.txt", "w");
	fprintf(fp_Spit, "This file spits out certain variables of interest\nSM Masud Rana\n--\n");
	fprintf(fp_Spit, "senM.exe ran with arguments: [1] %s\t [2] %s\n--\n", argv[1], argv[2]);
	fprintf(fp_Spit, "INP file: %s\n", inp_file_ascii);
	fprintf(fp_Spit, "Scada db config file: %s\n", db_config_file);
	fprintf(fp_Spit, "Model file: %s\n", xd_model_file);
	fprintf(fp_Spit, "Act Dem file: %s\n", actual_demand_file);
	fprintf(fp_Spit, "ClusterList file: %s\n", cluster_list_file);
	fprintf(fp_Spit, "is_dev_demand:\t\t %d\n", is_dev_demand);
	fprintf(fp_Spit, "EMMC:iterations:\t %d\n", n_iter);
	fprintf(fp_Spit, "EMMC:estimation end time:\t\t %s\n", est_end_time); //sm.
	fprintf(fp_Spit, "EMMC:max chain:\t\t %d\n", mcmc_chain_size);
	fprintf(fp_Spit, "EMMC:burn in:\t\t %d\n", burn_in);
	fprintf(fp_Spit, "EMMC:EstWinSize:\t\t %d\n", est_win_size); //sm.
	fprintf(fp_Spit, "EMMC:proposal_std:\t\t %f\n", prop_std); //sm.
	fprintf(fp_Spit, "EMMC:norm_std_ts:\t\t %f\n", norm_std_ts); //sm.
	fprintf(fp_Spit, "--\nModel file Parameters:\n----------------------\n"); //sm.
	fprintf(fp_Spit, "# of demand vec, dim =\t %d\n", dim); //sm.
	fprintf(fp_Spit, "# of seasons, ns     =\t %d\n", ns); //sm.
	fprintf(fp_Spit, "# of AR params, np     =\t %d\n", np);
	fprintf(fp_Spit, "mu[] = "); //sm.

	for (int i = 0; i < dim; ++i) {
		fscanf(xmf, "%lf", &mu[i]);
		fprintf(fp_Spit, "%.1f\t", mu[i]);
	}
	fprintf(fp_Spit, "\n"); //sm.

	int pattern_size;
	double *multiplier;
	double* real_demand; // real demands for a single time step

	real_demand = (double*)calloc(dim, sizeof(double));		// sm. real_demand[7]

															// tag: demand_multiplier
															//sm./ the following block of code is skipped if is_dev_deman = 0
	if (is_dev_demand) {
		fprintf(fp_Spit, "--\nSpitting out demand multipliers (tag:demand_multiplier):\n"); //sm.
		fscanf(xmf, "%d", &pattern_size);
		multiplier = (double*)calloc(pattern_size, sizeof(double));
		real_demand = (double*)calloc(dim, sizeof(double));		// sm. real_demand[7]
		for (int i = 0; i<pattern_size; ++i) {
			fscanf(xmf, "%lf", &multiplier[i]);
			fprintf(fp_Spit, "%f\n", multiplier[i]); //sm.
		}
	}

	//sm./ initial phi for phi11,phi12,phi13...
	for (int i = 0; i<dim*dim*np; ++i) {
		fscanf(xmf, " %lf", &phi0[i]);
	}

	fprintf(fp_Spit, "--\nInitial phi matrices (tag: initial phi).\n"); //sm.

																		//tag: initial phi
	int countsm = 0; // sm. counter 
	for (int j = 0; j < np; j++)
	{
		fprintf(fp_Spit, "Phi%d:\n", j);
		for (int i = 0; i < dim; i++) {
			for (int k = 0; k < dim; k++)
			{
				fprintf(fp_Spit, "%3.2f\t", phi0[countsm]); //sm.
				countsm++;
			}
			fprintf(fp_Spit, "\n");
		}
		fprintf(fp_Spit, "--\n");
	}
	fprintf(fp_Spit, "Successfully printed initial Phi Matrix (tag = initial phi)(count = %d, dim*dim*np = %d)\n--\n", countsm, dim*dim*np);

	//sm. Scanning in the initial covariance matrix
	//tag:cov initial
	for (int i = 0; i<dim*dim; ++i) {
		fscanf(xmf, " %lf", &cov0[i]);
	}

	countsm = 0;
	fprintf(fp_Spit, "Cov Matrix:\n");
	for (int i = 0; i < dim; i++) {
		for (int k = 0; k < dim; k++)
		{
			fprintf(fp_Spit, "%3.2f\t", cov0[countsm]); //sm.
			countsm++;
		}
		fprintf(fp_Spit, "\n");
	}
	fprintf(fp_Spit, "\n");
	fprintf(fp_Spit, "Successfully printed initial Cov Matrix (tag:cov initial)\n(count = %d, dim*dim = %d)\n--\n", countsm, dim*dim);

	// est_win_size = 168 form Net1 from model file
	double *xd0; //SM. initial water demand that we put in the MOD file
	xd0 = (double*)calloc(dim*est_win_size, sizeof(double));
	for (int i = 0; i<dim*est_win_size; ++i) {
		fscanf(xmf, "%lf", &xd0[i]);
	}

	countsm = 0;
	fprintf(fp_Spit, "Spitting out Initial Water Demand:\n");
	for (int i = 0; i <est_win_size; i++) {
		for (int k = 0; k < dim; k++)
		{
			fprintf(fp_Spit, "%.3f\t", xd0[countsm]); //sm.
			countsm++;
		}
		fprintf(fp_Spit, "\n");
	}
	fprintf(fp_Spit, "--\n");
	fprintf(fp_Spit, "Successfully printed initial demand values (count = %d, dim*est_win_size = %d)\n\n", countsm, dim*est_win_size);

	fclose(xmf);		// model file closed

						// row major to column major
	for (int ip = 0; ip<np; ++ip)
		for (int ir = 0; ir<dim; ++ir)
			for (int ic = 0; ic<dim; ++ic)
			{
				phi[ip*dim*dim + ic*dim + ir] = phi0[ip*dim*dim + ir*dim + ic];
			}
	//masud: phi is obtaining the transpose of the phi we provided in the MOD file.

	for (int ir = 0; ir<dim; ++ir)
		for (int ic = 0; ic<dim; ++ic) {
			cov[ic*dim + ir] = cov0[ir*dim + ic];
		}
	//masud: cov is obtaining the transpose of the cov we provided in the MOD file.

	VarModel_Err vmec = VarModel_new(dim, ns, s, p, &vm);
	VarModel_estMu(vm, mu, dim, 1);
	VarModel_setPhi(vm, phi, dim*dim*np);
	VarModel_Err vmec_cov = VarModel_setCov(vm, cov, dim*dim);

	//if (mode == 0) VarModel_dump(vm);
	VarModel_dump(vm);

	////----+++++++++++ACTUAL DEMAND FOR FIRST LOOKBACK HOURS++++++++++++++++++++++++++++++++++++++++++++++------------////////////////////
	//masud: 
	// First lookback hours data are required for demand likelihood calculations using a normal distribution (check std of normal dist below:)
	// Jinduan used to copy the last 26 hrs and used it as 
	
	double* dem_actual = (double*)calloc(vm->lb*dim, sizeof(double));
	FILE *fptemp; // "demand_actual.txt" -- read
	int fileOpenError = fopen_s(&fptemp, actual_demand_file, "r");
	if (fileOpenError)
	{
		printf("\n\nCannot open Actual Demand File (\"%s\")\n",actual_demand_file);
		return 1000;
	}
	for (int ii = 0; ii < vm->lb; ii++)
		for (int jj = 0; jj < dim; jj++)
			fscanf(fptemp, "%lf", &dem_actual[ii*dim + jj]);
	fclose(fptemp);

	//masud: temporary variable for demand likelihood calc for first 26hrs
	double* demands_ts = (double*)calloc(dim, sizeof(double));

	////----+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++------------////////////////////

	FILE *fp167;  //"debug_varm.txt"
	fopen_s(&fp167, "debug_varm.txt", "w");
	VarModel_dump2(vm, 0, &fp167);
	fclose(fp167);

	free(s); free(p);
	free(phi); free(phi0); free(cov0); free(cov);

	//masud:
	//cluster : (6-7-2016) commenting the following line for clustering
	if (vm->n != sv->_nXd) {
		printf("Dimensions of the demand model and the network/solver do not match. (Masud: OK if you are clustering!)\n");
		//return 110;
	}

	if (est_win_size < vm->lb) {
		printf("\n\nSize of demand data for initializing is shorter than demand model lookback.\n");
		printf("Only E-Step will be performed.\n");
		mode = 2;
		//return 120;
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
	Pop_Err pec = Pop_new(mcmc_chain_size, dim, 1, &pop);

	// Loading scada panel data
	//tag:scada data
	printf("[0]Loading scada panel data from the channels for estimation...\n");
	printf("  Time interval %d s, Current time (Year month day) in EM: %s. Using "
		"previous %d snapshots.\n",
		ds->dt, est_end_time, est_win_size);
	Tstamp tcur;  // estimation end time for EM, usually the current time stamp
	sscanf(est_end_time, "%d %d %d %d %d %d",
		&tcur.year, &tcur.month, &tcur.day,
		&tcur.hour, &tcur.minute, &tcur.second);
	tcur.fraction = 0;

	double* scada; //= size(scada) = 7*168
	int n_ch = ds->n_chan; //number of channels
	int ws = est_win_size; //window size = 168
	scada = (double*)calloc(n_ch * ws, sizeof(double));
	dsec = ds->fillSnapshots(tcur, ws, scada); // 168 data points per channel 

	fprintf(fp_Spit, "--\nSpitting out SCADA data input from the variable scada (tag:scada data):\n");
	countsm = 0;
	for (int j = 0; j < ws; j++) {
		for (int i = 0; i < n_ch; i++) {
			fprintf(fp_Spit, "%6.1f\t", scada[countsm++]);
		}
		fprintf(fp_Spit, "\n");
	}
	fprintf(fp_Spit, "--\n");

	fclose(fp_Spit);

	if (mode == 0) {
		printf("SenM checked out successfully.");
		goto END;
	}

	//demand data buffer for var model look-back and estimation
	int buf_size = vm->lb + ws; // look back + window of est = 168+26 = 194
	double* buf = (double*)calloc(buf_size*dim, sizeof(double)); //buf[194][7] 194*7 = 1358 == buffer contains initial data xd0

																 //load init demand - for a continuously running program this may
																 // be the previous demand estimates
																 //masud: copying the last 26 hours (latest) of initial demand from xd0 to buf
	memcpy(buf, &xd0[(est_win_size - vm->lb)*dim],
		sizeof(double)*vm->lb*dim); // lb*dim = 26*7 = 182
									//sm. copying the entire set of initial demand from xd0 to buf
	memcpy(&buf[vm->lb*dim], xd0, sizeof(double)*ws*dim); //1176


	FILE *fp163;	//Masud: "Debug_Likelihood.txt" -- //file 2016 -2; for debug outputs
	fopen_s(&fp163, "Debug_Likelihood.txt", "w");
	fprintf(fp163, "\nHyd L       Demand L   Total L     Diff in Mu and Phi\n");

	FILE *fp164;	//Masud: "Debug_FlowCalculations.csv"
	fopen_s(&fp164, "Debug_FlowCalculations_Avg.csv", "w");
	FILE *fpFlowChain1;	//Masud: "Debug_FlowCalculations.csv"
	fopen_s(&fpFlowChain1, "Debug_FlowCalc_t49.csv", "w");
	FILE *fpFlowChain2;	//Masud: "Debug_FlowCalculations.csv"
	fopen_s(&fpFlowChain2, "Debug_FlowCalc_t107.csv", "w");
	FILE *fpFlowChain3;	//Masud: "Debug_FlowCalculations.csv"
	fopen_s(&fpFlowChain3, "Debug_FlowCalc_t137.csv", "w");

	FILE *fp_LLchain1;
	fopen_s(&fp_LLchain1, "Debug_FlowChain_t49.csv", "w");
	FILE *fp_LLchain2;
	fopen_s(&fp_LLchain2, "Debug_FlowChain_t107.csv", "w");
	FILE *fp_LLchain3;
	fopen_s(&fp_LLchain3, "Debug_FlowChain_t137.csv", "w");

	FILE *fp_buffer; //"Debug_Buffer2.txt"
	fopen_s(&fp_buffer, "Debug_Buffer2.txt", "w");

	fopen_s(&fp167, "Debug_varm.txt", "a");

	////----+++++++++++++++++++++++ DEMAND GROUPING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---||||||||||
	// MASUD -- 6-8-2016

	// todo 1: implement as user choice in the model file

	//is_cluster = 1; //masud: flag to indicate if clustering of demand nodes will be employed or not. 
						//If true, initial guess should be multiplier rather than actual demand

	int* cluster_id = (int*)calloc(net->Njuncs, sizeof(int));
	FILE *fp_cluster;
	int ferr = fopen_s(&fp_cluster, cluster_list_file, "r");
	if (ferr)
	{
		printf("\n\nCannot open \"%s\"\n\n",cluster_list_file);
	}

	for (int i = 0; i < net->Njuncs; i++)
		int err = fscanf(fp_cluster, "%d", &cluster_id[i]);
	fclose(fp_cluster);

	//FILE *fp_llchain; // LikelihoodChain.txt
	//fopen_s(&fp_llchain, "LikelihoodChain.txt", "w");
	//fprintf(fp_llchain, "\nHyd L       Demand L   Total L\n");

	////----+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++------------////////////////////

	//init random number generator
	Rvgs* rng = new Rvgs(rng_seed); // random number generator
	printf("\nStarting EM-Cycle\n");
	bool print_all_PQ = false;

	clock_t time_beginning = clock();
//--------------------------------// EM cycle //-------------------------------------------------------------------------------------------//
	for (int iter = 0; iter < n_iter; ++iter)
	{
		printf("\nEM iteration %d (Percent acceptance of MC sample shown below): \n", iter);
		double t_pct_acc = 0; //sm. total percent accepted

							  // E-step: t = 0:167
		for (int t = 0; t < ws; ++t)
		{
			//#define POPH(pop, iw, im, id)   (pop->h[iw][(id) + (im)*pop->d1*pop->d2])
			//first sample
			Pop_reset(pop);
			//memcpy(&POPH(pop, 0, 0, 0), mu , sizeof(double)*dim);
			memcpy(&POPH(pop, 0, 0, 0), &buf[(vm->lb + t)*dim], //masud: h[][0] has only one hour of data... starting from t=0...through ws
				sizeof(double)*dim);

			// Masud: pop->h[][] has buf[(26+t)*7] to buf[(26+t)*7 + 7] 
			// == basically pop-->h[][] has one hour worth of initial data starting from the most 
			// current time minus lookback length
			// as t increases in this loop (up to ws)

			// MCMC
			double cur_ll = 0;  //current (and Hydraulic) log-likelihood
			double hyd_likelihood = 0; // hydraulic likelihood to be calculated from hydraulics and observed scada data
			double dem_likelihood = 0; // demand likelihood to be calculated from Time Series or Normal distribution for the first lookback hours

			double last_ll = -DBL_MAX;  //ll from the previous sample
			int n_acc = 1;
			int n_rej = 0;
			////////////////////////////////////// BEGINNING OF E-STEP //////////////////////////////////////////


			for (int im = 0; im < mcmc_chain_size; ++im) //im = 1000 chain
			{
				if (im == burn_in)
				{
					n_acc = 1;
					n_rej = 0;
				}

				double* cur_xd = &POPH(pop, 0, im, 0); // masud: cur_xd[7] = h[][im*dim] : h[][] has initial demand from MOD
													   // masud: from h[][im*dim] to h[][im*dim +7] there is only one hour worth of data... rest is zeros.

													   // hydraulic likelihood
				if (is_dev_demand)
				{
					sv->run(real_demand, dim, &scada[t*n_ch], n_ch);
				}
				else if (is_log_demand)
				{
					sv->runlogd(cur_xd, dim, &scada[t*n_ch], n_ch);
				}
				else
				{
					//masud: this is where the hydraulic solver is used:
					//sv->run(cur_xd, dim, &scada[t*n_ch], n_ch);
					int dim_cluster = sv->_nXd; // Masud: lying to run2() by saying demand variables are more than what they actually are.
					sv->run2(cur_xd, dim_cluster, &scada[t*n_ch], n_ch, cluster_id, is_cluster);
				}

				bool smbool = false;
				//if (t == 6 && im == 0)
				//smbool = true;

				//************************* PRINTING FLOW CHAIN AT PARTICULAR HOURS ****************************//
				// MASUD: WARNING: Hard coded hours. Remove the hard coded hours to replace with user-defined hours to prevent crash.
				//sv->logL(&scada[t*n_ch], n_ch, &hyd_likelihood, false, smbool);
				 //print pressures and flow
				if (t == 49)
				{
					print_all_PQ = true;
					sv->logL2(&scada[t*n_ch], n_ch, &hyd_likelihood, &fpFlowChain1, print_all_PQ);
				}
				else if (t == 107)
				{
					print_all_PQ = true;
					sv->logL2(&scada[t*n_ch], n_ch, &hyd_likelihood, &fpFlowChain2, print_all_PQ);
				}
				else if (t == 137)
				{
					print_all_PQ = true;
					sv->logL2(&scada[t*n_ch], n_ch, &hyd_likelihood, &fpFlowChain3, print_all_PQ);
				}
				else
				{
					print_all_PQ = false;
					sv->logL2(&scada[t*n_ch], n_ch, &hyd_likelihood, &fp164, print_all_PQ);
				}
					
				
				
				//************************* PRINTING FLOW CHAIN END ****************************//

				//masud:logL2 gets hydraulic likelyhood details are a particular time for a particular iteration
				/*if(t==130)
				sv->logL2(&scada[t*n_ch], n_ch, &hyd_likelihood, &fp164);
				else
				sv->logL(&scada[t*n_ch], n_ch, &hyd_likelihood, false, smbool);*/

				// compute demand likelihood
				// masud: buf variable is updated here with current sample
				//double dm_ll=0;
				memcpy(&buf[(vm->lb + t)*dim], cur_xd, sizeof(double)*dim);

				//%%%%%%%%%%%%%%%%% Time Serie Likelihood Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				if (t >= vm->lb)
				{
					VarModel_logLikelihood(vm, &buf[t*dim], dim, vm->lb + 1, &dem_likelihood);
				}
				else
				{
					for (int ii = 0; ii < dim; ii++)
						demands_ts[ii] = buf[vm->lb*dim + t*dim + ii];

					double temp_ts_likelihood = 0;
					double ts_likelihood = 0; //time series likelihood
					//double ts_sigma = .8; // this sigma is only to calculate TS prior for the first 26 hrs
					//norm_std_ts
					for (int ii = 0; ii < dim; ii++)
					{
						//temp_ts_likelihood = pdfNorm(demands_ts[ii], dem_actual[t][ii], norm_std_ts);
						temp_ts_likelihood = pdfNorm(demands_ts[ii], dem_actual[t*dim + ii], norm_std_ts);
						ts_likelihood = ts_likelihood + log(temp_ts_likelihood);
					}

					dem_likelihood = ts_likelihood;
				}
				//%%%%%%%%%%%%%%%%% TS Likelihood END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				cur_ll = dem_likelihood + hyd_likelihood;
				
				if (t == 49)
				{
					fprintf(fp_LLchain1, "%f,%f,%f\n", dem_likelihood, hyd_likelihood, cur_ll);
				}
				else if (t == 107)
				{
					fprintf(fp_LLchain2, "%f,%f,%f\n", dem_likelihood, hyd_likelihood, cur_ll);
				}
				else if (t == 137)
				{
					fprintf(fp_LLchain3, "%f,%f,%f\n", dem_likelihood, hyd_likelihood, cur_ll);
				}

				
				/////////////////// ACCEPT-REJECT using Metropolis-Hastings algorithm) ///////////////////////
				int poorer_likelihood = cur_ll < last_ll; // TRUE if current hyd likelyhood is worse
				int rejection_prob = rng->Random() > exp(cur_ll - last_ll);		// ratio of current over previous likelihoods

				//MH accept-reject // tag: MH_alg
				//masud: Replacement of the original line above for debugging.
				if (poorer_likelihood && rejection_prob)	//masud: poorer_likelihood is redundant?? -- YES since Random() > exp(cur_ll - last_ll) is always false 
				{																						// for cur_ll > last_ll				
					// reject, roll-back to the previous sample
					memcpy(cur_xd, &POPH(pop, 0, im - 1, 0), sizeof(double)*dim);

					cur_ll = last_ll;
					++n_rej;
					//printf("x");
				}
				else {
					++n_acc;
					//printf("-");
				}

				//advance chain pointer 
				pop->p[0]++;
				last_ll = cur_ll;

				//masud: this is where markov chan is created by adding a normally distributed random noise
				if (im == mcmc_chain_size - 1) break;
				else
				{//propose a new demand panel
					for (int id = 0; id<dim; ++id)
					{
						POPH(pop, 0, im + 1, id) = POPH(pop, 0, im, id) + rng->Normal(0, prop_std);

						//printf("NoChain ");
						//POPH(pop, 0, im + 1, id) = POPH(pop, 0, im, id);
					}
				}

				/*if (im % 1000 == 0)
				{
				double pct_acc = 100.0 * n_acc / (n_acc + n_rej);
				printf("(%3.0f%%) ", pct_acc);
				}*/

			}//End of for (int im = 0; im < mcmc_chain_size ; ++im) //im = 10000 chain
			 //End of for (int im = 0; im < mcmc_chain_size ; ++im) //im = 10000 chain
			 //End of for (int im = 0; im < mcmc_chain_size ; ++im) //im = 10000 chain
			 //mcmc for one hour ended
			 ///////////////////// End of MC Chain ///////////////////////////////////////////////////////////////////////
			
			double pct_acc = 100.0 * n_acc / (n_acc + n_rej);
			printf("%3.0f%% ", pct_acc);

			//Pop_calc(pop);
			//Pop_report(pop);
			//Pop_writeout(pop, estep_outfile, NULL, mcmc_chain_size, t);

			//if(iter == 4 || iter == 9)
			Pop_writeout2(pop, estep_outfile, NULL, mcmc_chain_size, t, iter);

			Pop_Err perr = Pop_mean(pop, &buf[(vm->lb + t)*dim], dim, burn_in); //update demand estimates with sample means
			t_pct_acc += pct_acc;

			if (perr) return 301;
		} // End of loop: for (int t = 0; t < ws; ++t)
		  // End of loop: for (int t = 0; t < ws; ++t)
		  // End of loop: for (int t = 0; t < ws; ++t)

		t_pct_acc /= ws;
		printf("\nAverage Acceptance of all Hours (pro std = %.4f): %f%%\n\n", prop_std, t_pct_acc);
		////////////////////////SM.-----------END OF E-Step---------------------------------------------------///////

		////////////////////////SM.-----------M -Step ---- M -Step -----M -Step ------M -Step  ---------------///////
		double diff; // masud: difference of  mu
		double mdll; //TS model (VAR) likelihood
		//%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//masud: debug 4
		//M-step: re-estimate demand model
		//VarModel_Err vmerr = VarModel_estimate(vm, buf, dim, buf_size, &diff, &mdll);

		if (mode != 2)	// mode =2 for E-step only
		{
			VarModel_Err vmerr = VarModel_estimate2(vm, buf, dim, buf_size, &diff, &mdll);
			if (vmerr) return 401;
			printf("\nVAR Model updated in M step (Error code: %d). Iteration %d\n", vmerr, iter);
			printf("\n");
			//VarModel_dump(vm);
			VarModel_dump2(vm, iter + 1, &fp167); //hb just prints the cov matrix
		}

		//masud: for log demand
		if (is_log_demand) {
			printf("exp(mu):\n");
			for (int i = 0; i<dim; ++i) {
				printf("%6.3f ", exp(vm->mu[i]));
			}
		}
		double mhll = 0; // mean hydraulic likelihood

		//sm./ tag:bug -- the program tries to use the var multiplier but it hasn't been initialized.
		// the following block produces the numbers in the console... if the format %d, %f, %f...
		for (int t = 0; t < ws; ++t)
		{
			double mhllt = 0; //masud: temporary hydraulic likelyhood
							  /*{
							  for (int i=0; i<dim; ++i)
							  real_demand[i] = multiplier[t%pattern_size]*mu[i]
							  + buf[(vm->lb+t)*dim+i];
							  }*/
			for (int i = 0; i<dim; ++i)
			{ //masud: this block is code is modified to get rid of the bug (5/12/2016) see above for the original block
				if (is_dev_demand)
				{
					real_demand[i] = multiplier[t%pattern_size] * mu[i]
						+ buf[(vm->lb + t)*dim + i];
				}
				else
				{
					real_demand[i] = buf[(vm->lb + t)*dim + i];
				}
			}

			//masud : sv solves the network hydraulics with the mean estimates of demand
			//sv->run(real_demand, dim, &scada[t*n_ch], n_ch);
			int dim_cluster = sv->_nXd;
			sv->run2(real_demand, dim_cluster, &scada[t*n_ch], n_ch, cluster_id, is_cluster);

			bool print_PQ = true; //masud: printing observed and calculated pressure/flow and likelihood for mean demand
			sv->logL2(&scada[t*n_ch], n_ch, &mhllt, &fp164, print_PQ);	//masud: logL 
																		//sv->logL(&scada[t*n_ch], n_ch, &mhllt, true);	//masud: logL produces the %d %f %f output to console. 
			mhll += mhllt;
		}

		printf("Sum of changes in parameters:\n "
			"Sqaured diff in mu and phi %.2f, VAR likelihood: %.2f, Hyd Likelihood: %.2f, total likelihood:%.2f\n",
			diff, mdll, mhll, mdll + mhll);
		fprintf(fp163, "%f %f %f %f\n", mhll, mdll, mdll + mhll, diff);

	}
	// ////////////////////////SM.-----------END OF M-Step---------------------------------------------------///////
	//sm. End of : for (int iter = 0; iter < n_iter; ++iter)
	//sm. End of : for (int iter = 0; iter < n_iter; ++iter)
	//sm. End of : for (int iter = 0; iter < n_iter; ++iter)

	clock_t time_end = clock();
	double time_req = (double)(time_end - time_beginning) / CLOCKS_PER_SEC;
	printf("\nTotal Time required for %d EM Cycles is %.1f sec\nAverage time per cycle = %.1f sec\n\n", n_iter, time_req, time_req / n_iter);
	fprintf(fp167,"\n\n\nTotal Time required for %d EM Cycles is %.1f sec\nAverage time per cycle = %.1f sec\n\n", n_iter, time_req, time_req / n_iter);

SM_TERMINATE:
	printf("\nSuccessfully Finished E-M Iterations\n\nCleaning up memory...");
	//masud: commenting the following:
	//print_mat(buf, dim, buf_size, "Buffer");	//printing buf
	//fclose(fp_llchain);
	fclose(fp_buffer);
	fclose(fp163);
	fclose(fp164);
	fclose(fpFlowChain1);
	fclose(fpFlowChain2);
	fclose(fpFlowChain3);

	fclose(fp_LLchain1);
	fclose(fp_LLchain2);
	fclose(fp_LLchain3);
	//fclose(fp165);
	fclose(fp167);
	// masud: for early termination of simulation: 

	//----------Masud: freeing memory of new variables
	free(demands_ts);
	free(dem_actual);
	free(cluster_id);
	//--------------------------

	free(buf);
	free(real_demand);
	free(xd0);
	free(mu);
	free(scada);
	delete rng;

END:
	// clean up

	VarModel_del(&vm);
	Pop_del(&pop);


	// Following Crushes
	//delete sv;		//masud: delete sv causes access violation error
	//delete ds; // masud: this causes an access violation error

	// TODO: Network destructor
	//delete net;

	iniparser_freedict(dict);
	printf("SUCCESS!!\n\n");

	return 0;
}







