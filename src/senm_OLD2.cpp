// EMMC program of senm system
// jinduan:
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include <iostream>

#include "iniparser.h"
#include "rvgs.h"

#include "SenMCoreIncs.h"

#include "Network.h"
#include "DataSource.h"
#include "Solver.h"
#include "Population.h"
#include "VarModel.h"

using namespace std;

int main(int argc, char** argv) {
	// argv[1] should be "-estimate" for mode 1
	// argv[2] should be the config file i.e., senm.ini
	FILE *fp16; // fp 2016
	

	//sm. debug file: I am going to output different variables in a text file after reading
	fopen_s(&fp16, "spitOut.txt", "w");
	
	fprintf(fp16, "This file spits out certain variables of interest\nSM Masud Rana\n--\n");

	USES_CONVERSION;

	if (argc != 3) {cout<< "Missing Input Parameters. Rerun the program with two parameters.\n"<<endl; return 1; }
	
	printf("argc = %d \n argv = %s, %s\n--\n", argc, argv[1], argv[2]);

	int mode = 0; // data integrity check
	if (strcmp(argv[1], "-estimate") == 0) {
		mode = 1; // emmc
	}
	//smmr: strcmp() returns 0 when strings exactly match!

	printf("[0]Loading config file %s ... \n", argv[2]);
	dictionary* dict;
	dict = iniparser_load(argv[2]);
	if (dict == NULL) {
		return 2;
	}

	fprintf(fp16, "senM.exe ran with arguments: [1] %s\t [2] %s\n--\n",argv[1],argv[2]);

	char* inp_file_ascii = iniparser_getstring(dict, "Net:inp_file", "");
	if (strcmp(inp_file_ascii, "") ==0) {
		return 3;
	}
	
	fprintf(fp16, "INP file: %s\n", inp_file_ascii);
	
	TCHAR inp_file_tchar[MAX_FILE_NAME_SIZE];
	_tcscpy(inp_file_tchar, A2T(inp_file_ascii));

	char* db_config_file = iniparser_getstring(dict, "Scada:db_conn_file", "");
	if (strcmp(db_config_file, "") ==0) {
		return 5;
	}
	fprintf(fp16, "Scada db config file: %s\n", db_config_file);

	char* xd_model_file = iniparser_getstring(dict, "XDemands:model_file", "");
	if (strcmp(xd_model_file, "") ==0) {
		return 5;
	}
	fprintf(fp16, "Model file: %s\n", xd_model_file);

	

	int is_log_demand = iniparser_getint(dict, "XDemands:is_log_demand", 0);
	int is_dev_demand = iniparser_getint(dict, "XDemands:is_dev_demand", 0);

	fprintf(fp16, "is_dev_demand:\t\t %d\n", is_dev_demand);

	int n_iter = iniparser_getint(dict, "EMMC:em_iterations", 0);
	if (n_iter <= 0) {
		return 106;
	}
	fprintf(fp16, "EMMC iterations:\t %d\n", n_iter);

	int mcmc_chain_size = iniparser_getint(dict, "EMMC:max_chain_size", 0);
	if (mcmc_chain_size <= 0) {
		return 101;
	}
	
	fprintf(fp16, "EMMC:max chain:\t\t %d\n", mcmc_chain_size);

	int burn_in = iniparser_getint(dict, "EMMC:burn_in", 0);
	if (burn_in <= 0 || burn_in > mcmc_chain_size) {
		return 107;
	}
	fprintf(fp16, "EMMC burn in:\t\t %d\n", burn_in);
/*
	int init_xd_size = iniparser_getint(dict, 
		"EMMC:initial_demand_data_length", 0);
	if (init_xd_size <= 0) {
		return 102;
	}
	*/

	int est_win_size = iniparser_getint(dict, 
		"EMMC:estimation_win_size", 0); // est_win_size = 168 for Net1
	if (est_win_size <= 0) {
		return 104;
	}
	fprintf(fp16, "Est window size:\t\t %d\n", est_win_size); //sm.

	double prop_std = iniparser_getdouble(dict, 
		"EMMC:proposal_std", 0);
	if (prop_std <=0) {
		return 103;
	}
	fprintf(fp16, "EMMC:proposal_std:\t\t %f\n", prop_std); //sm.

	char* est_end_time = iniparser_getstring(dict, 
		"EMMC:estimation_end_time", "");
	if (strcmp(est_end_time, "") == 0) {
		return 13;
	}
	fprintf(fp16, "EMMC:estimation end time:\t\t %s\n", est_end_time); //sm.

	char* estep_outfile = iniparser_getstring(dict, "EMMC:estep_outfile", "");
	if (strcmp(estep_outfile, "") == 0) return 14;

	char* os_outfile = iniparser_getstring(dict, "EMMC:os_debug_outfile", "");
	if (strcmp(os_outfile, "") == 0) os_outfile = NULL;

	int output_interval = 80;
//		iniparser_getint(dict, "EMMC:data_output_interval", 0);
//	if (output_interval == 0) return 15;


	
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

	//tag:database connection
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
	int dim; // sm. number of demand variables
	int ns;	// sm. number of seasons
	int* s;
	int* p;
	int np;	// sm. total number of AR params
	double *mu;
	double *phi0, *phi; // row- and column- major
	double *cov0, *cov;
	FILE* xmf;

	//tag:modelFile start
	xmf = fopen(xd_model_file, "r");
	fscanf(xmf, "%d", &dim);	//sm./ dim == first line of the .mod file. = number of demand var
	fscanf(xmf, "%d", &ns);		//sm./ dim == second line of the .mod file. = number of seasons

	
	fprintf(fp16, "--\nModel file Parameters:\n----------------------\n"); //sm.
	fprintf(fp16, "(tag:modelFile start) # of demand vec, dim =\t %d\n", dim); //sm.
	fprintf(fp16, "(tag:modelFile start) # of seasons, ns     =\t %d\n", ns); //sm.

	if (dim<=0 || ns <=0) return 8;

	s = (int*)malloc(sizeof(int)*(ns+1));	// ns = 1 for Net1
	p = (int*)malloc(sizeof(int)*(ns+1));

	//sm./ ns = 1 for NET1
	for (int is=0; is<=ns; ++is) {
		fscanf(xmf, "%d", &s[is]);
	}

	np = 0;
	for (int is=0; is<=ns; ++is) {
		fscanf(xmf, "%d", &p[is]);
		np += p[is];
	}//sm./ total number of AR parameters

	fprintf(fp16, "(Line 206, senm.cpp) # of seasons, np     =\t %d\n", np);

	mu = (double*)calloc(dim, sizeof(double)); //sm./ dim = 7 for Net1
	cov0 = (double*)calloc(dim*dim, sizeof(double));
	cov = (double*)calloc(dim*dim, sizeof(double));
	phi0 = (double*)calloc(dim*dim*np, sizeof(double));
	phi = (double*)calloc(dim*dim*np, sizeof(double));

	//tag : mu
	fprintf(fp16, "(tag: mu, senm.cpp) mu[] = "); //sm.

	for (int i = 0; i < dim; ++i ) {
		fscanf(xmf, "%lf", &mu[i]);
		fprintf(fp16, "%.1f\t", mu[i]);
	}
	fprintf(fp16, "\n"); //sm.

	
	
	int pattern_size;
	double *multiplier;
	double* real_demand; // real demands for a single time step
	
	real_demand = (double*)calloc(dim, sizeof(double));		// sm. real_demand[7]

	// tag:demand multiplier
	//sm./ the following block of code is skipped if is_dev_deman = 0
	if (is_dev_demand) {
		fprintf(fp16, "--\nSpitting out demand multipliers (tag:demand multiplier):\n"); //sm.
		fscanf(xmf, "%d", &pattern_size);
		multiplier = (double*)calloc(pattern_size, sizeof(double));
		real_demand = (double*)calloc(dim, sizeof(double));		// sm. real_demand[7]
		for (int i=0; i<pattern_size; ++i) {
			fscanf(xmf, "%lf", &multiplier[i]);
			fprintf(fp16, "%f\n", multiplier[i]); //sm.
		}
	}


	//sm./ initial phi for phi1,phi2,phi3...
	for (int i=0; i<dim*dim*np; ++i) {
		fscanf(xmf, " %lf", &phi0[i]);
	}

	fprintf(fp16, "--\nInitial phi matrices (tag: initial phi).\n"); //sm.

	//tag: initial phi
	int countsm = 0; // sm. counter 
	for (int j = 0; j < np; j++)
	{
		fprintf(fp16, "Phi%d:\n", j);
		for (int i = 0; i < dim; i++) {
			for (int k = 0; k < dim; k++)
			{
				fprintf(fp16, "%3.2f\t", phi0[countsm]); //sm.
				countsm++;
			}
			fprintf(fp16, "\n");
		}
		fprintf(fp16, "--\n");
	}
	fprintf(fp16, "Successfully printed initial Phi Matrix (tag = initial phi)(count = %d, dim*dim*np = %d)\n--\n",countsm,dim*dim*np);

	//sm. Scanning in the initial covariance matrix
	//tag:cov initial
	for (int i =0; i<dim*dim; ++i) {
		fscanf(xmf, " %lf", &cov0[i]);
	}

	countsm = 0;
	fprintf(fp16, "Cov Matrix:\n");
	for (int i = 0; i < dim; i++) {
		for (int k = 0; k < dim; k++)
		{
			fprintf(fp16, "%3.2f\t", cov0[countsm]); //sm.
			countsm++;
		}
		fprintf(fp16, "\n");
	}
	fprintf(fp16, "\n");
	fprintf(fp16, "Successfully printed initial Cov Matrix (tag:cov initial)\n(count = %d, dim*dim = %d)\n--\n", countsm, dim*dim);
	
	// est_win_size = 168 form Net1 from model file
	double *xd0; //SM. initial water demand that we put in the MOD file
	xd0 = (double*)calloc(dim*est_win_size, sizeof(double));
	for (int i=0; i<dim*est_win_size; ++i) {
		fscanf(xmf, "%lf", &xd0[i]);
	}

	/*/--------Debug----------tag:init demand----------------------------//
	// SM. I am using the following for loop for debugging only to override the initial demand reading 
	// from the MOD file. 
	double init_demand_sm = 30.5;
	for (int i = 0; i < dim*est_win_size; i++) {
		xd0[i] = init_demand_sm;
		//xd0[i] = xd0[i]/30;
	}
	//--------Debug--------------------------------------*/

	countsm = 0;
	fprintf(fp16, "Spitting out Initial Water Demand:\n");
	for (int i = 0; i <est_win_size; i++) {
		for (int k = 0; k < dim; k++)
		{
			fprintf(fp16, "%4.1f\t", xd0[countsm]); //sm.
			countsm++;
		}
		fprintf(fp16, "\n");
	}
	fprintf(fp16, "--\n");
	fprintf(fp16, "Successfully printed initial demand values (count = %d, dim*est_win_size = %d)\n\n", countsm, dim*est_win_size);
	
	fclose(xmf);		// model file closed

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
	VarModel_estMu(vm, mu, dim, 1);
	VarModel_setPhi(vm, phi, dim*dim*np);
	VarModel_Err vmec_cov = VarModel_setCov(vm, cov, dim*dim);

	//if (mode == 0) VarModel_dump(vm);
	VarModel_dump(vm);

	free(s); free(p);
	free(phi); free(phi0); free(cov0); free(cov);

	if (vm->n != sv->_nXd) {
		printf("Dimensions of the demand model and the network/solver do not match.\n");
		return 110;
	}

	if (est_win_size < vm->lb ) {
		printf("size of demand data for initializing is shorter than demand model lookback.\n");
		return 120;
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

	/*
	// Masud: following block moved after scada data reading.
	if (mode == 0) {
		printf("SenM checked out successfully.");
		goto END;
	}
	*/
	

	// Loading scada panel data
	//tag:scada data
	printf("[0]Loading scada panel data from the channels for estimation...\n");
	printf("  Time interval %d s, Current time (Year month day) in EM: %s. Using "
		   "previous %d snapshots.\n", 
		ds->dt,  est_end_time, est_win_size);
	Tstamp tcur;  // estimation end time for EM, usually the current time stamp
	sscanf(est_end_time, "%d %d %d %d %d %d", 
		&tcur.year, &tcur.month, &tcur.day, 
		&tcur.hour, &tcur.minute, &tcur.second);
	tcur.fraction = 0;
	
	double* scada; //= size(scada) = 7*168
	int n_ch = ds->n_chan; //number of channels
	int ws = est_win_size; //window size
	scada = (double*)calloc(n_ch * ws, sizeof(double));
	dsec = ds->fillSnapshots(tcur, ws, scada); // 168 data points per channel 

	fprintf(fp16, "--\nSpitting out SCADA data input from the variable scada (tag:scada data):\n");
	
	countsm = 0;
	for (int j = 0; j < ws; j++) {
		for (int i = 0; i < n_ch; i++) {
			fprintf(fp16, "%6.1f\t", scada[countsm++]);
		}
		fprintf(fp16, "\n");
	}
	fprintf(fp16, "--\n");

	fclose(fp16);

	if (mode == 0) {
		printf("SenM checked out successfully.");
		goto END;
	}

	//demand data buffer for var model look-back and estimation
	int buf_size = vm->lb + ws; // look back + window of est = 168+26 = 194
	double* buf = (double*)calloc(buf_size*dim, sizeof(double)); //buf[194][7] 194*7 = 1358 == buffer contains initial data xd0

	//load init demand - for a continuously running program this may
	// be the previous demand estimates
	
	//masud: copying the last 26 days (latest) of initial demand from xd0 to buf
	memcpy(buf, &xd0[(est_win_size - vm->lb)*dim], 
				sizeof(double)*vm->lb*dim); // lb*dim = 26*7 = 182
	
	/*printf("Debug: printing buffer:\n");
	for (int i = 0; i < buf_size; i++)
		printf("%5f  ", buf[i]);*/
	
	//sm. copying the entire set of initial demand from xd0 to buf
	memcpy(&buf[vm->lb*dim], xd0, sizeof(double)*ws*dim); //1176

	////SM. priting out buf:
	//printf("Debug: printing buffer:\n");
	//for (int i = 0; i < buf_size*dim; i++)
	//	printf("%5.1f  ", buf[i]);


	//init random number generator
	Rvgs* rng = new Rvgs(2); // random number generator

	FILE *fp162;	//Masud: file 2016 -2; for debug outputs
	fopen_s(&fp162, "Debug.txt", "w");
	
	double t_pct_acc = 0; //sm. total percent accepted
	// EM cycle
	for (int iter = 0; iter < n_iter; ++iter) 
	{
		printf("last average prct accepted: %f\n", t_pct_acc);
		printf("\nEM iteration %d (Percent acceptance of MC sample shown below. Current proposal standard dev = %4.1f): \n", iter,prop_std);
		
		t_pct_acc = 0; //sm. total percent accepted

		// E-step: sampling demands
		for (int t = 0; t < ws; ++t) 
		{
			//#define POPH(pop, iw, im, id)   (pop->h[iw][(id) + (im)*pop->d1*pop->d2])
			//first sample
			Pop_reset(pop);
			//memcpy(&POPH(pop, 0, 0, 0), mu , sizeof(double)*dim);
			memcpy(&POPH(pop, 0, 0, 0), &buf[(vm->lb+t)*dim], //masud: h[][0] has only one hour of data... starting from t=0...through ws
				sizeof(double)*dim);  
			
			// Masud: pop->h[][] has buf[(26+t)*7] to buf[(26+t)*7 + 7] 
			// == basically pop-->h[][] has one hour worth of initial data starting from the most 
			// current time minus lookback length
			// as t increases in this loop (up to ws)

			// MCMC
			double cur_ll = 0;  //current log-likelihood
			double last_ll = -DBL_MAX;  //ll from the previous sample
			int n_acc = 1;
			int n_rej = 0;
			
			for (int im = 0; im < mcmc_chain_size ; ++im) //im = 1000 chain
			{

				double* cur_xd = &POPH(pop, 0, im, 0); // masud: cur_xd[7] = h[][im*dim] : h[][] has initial demand from MOD
				int has_wrong_demand = 0;
				
				// masud: from h[][im*dim] to h[][im*dim +7] there is only one hour worth of data... rest is zeros.
				if (is_dev_demand) 
				{
					//tag:real demand	tag:here
					if (t < 50 && im<100) fprintf(fp162, "(tag:real demand)Printing real_demand[] || cur_xd[] mcm[i]=%d t=%d of %d(=ws)\n",im, t, ws);
					for (int i=0; i<dim; ++i) 
					{
						real_demand[i] = multiplier[t%pattern_size]*mu[i] + cur_xd[i];
						//sm. debug:
						if (t < 50 && im<100) {
							//fprintf(fp162, "%2.1f*%3.0f+%3.0f = %3.0f   ", multiplier[t%pattern_size], mu[i], cur_xd[i], real_demand[i]);
							fprintf(fp162, "%4.1f-%4.1f\t", real_demand[i],cur_xd[i]);
						}
					}
					//sm. debug:
					if(t<50 && im<100)
						fprintf(fp162, "\n");
				}
				
				//masud: tag:dom Q: why is he using initial guess's deviation from mu*pattern as a criteria for rejecting demand?
				//masud: why 0.2 is used as a multiplier? Current demand should not be that low!
				if (is_dev_demand) {
					for (int id = 0; id < dim; ++id) { //reject problematic demands
						double tp; //sm: temporary?
						tp = abs(cur_xd[id]) - multiplier[t%pattern_size]*mu[id] * 0.2;
						//tp = abs(real_demand[id]
						if (tp > 0) {
							//if ((rng->Random())/2 < tp) {
								has_wrong_demand = 1;
								break;
							//}
						}
					}
				} else {
					for (int id = 0; id < dim; ++id) { //reject problematic demands
						if (cur_xd[id] < xd0[t*dim+id] * 0.9 || 
							cur_xd[id] > xd0[t*dim+id] * 1.1) {
							has_wrong_demand = 1;
							break;
						}
					}
				}

				if (! has_wrong_demand) 
				{
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
						//dom: this also strengthens the hypothesis that deviation is used
						//masud: this is where the hydraulic solver is used:
						sv->run(cur_xd, dim, &scada[t*n_ch], n_ch);
					}
					sv->logL(&scada[t*n_ch], n_ch, &cur_ll);

					

					//compute demand likelihood
					// masud: buf variable is updated here for accepted samples
					double dm_ll=0;
					memcpy(&buf[(vm->lb+t)*dim], cur_xd, sizeof(double)*dim);
					VarModel_logLikelihood(vm, &buf[t*dim], dim, vm->lb+1, &dm_ll);

					cur_ll += dm_ll;
				}

				//MH accept-reject
				if (has_wrong_demand || cur_ll < last_ll && rng->Random() > exp(cur_ll - last_ll)) {
					// reject, roll-back to the previous sample
					memcpy(cur_xd, &POPH(pop, 0, im-1, 0), sizeof(double)*dim);
					
					cur_ll = last_ll;
					++n_rej;
					//printf("x");
				} else {
					++n_acc;
					//printf("-");
				}

				//advance chain pointer 
				pop->p[0]++;

				//if ((im+1) % output_interval == 0) printf("\n");	

				last_ll = cur_ll;

				//masud: this is where markov chan is created by adding a normally distributed random noise
				if (im == mcmc_chain_size - 1) break;
				else 
				{//propose a new demand panel
					for (int id=0; id<dim; ++id) 
					{
						POPH(pop, 0, im+1, id) = POPH(pop, 0, im, id) 
							+ rng->Normal(0, prop_std);
					}
				}

			} 
			//mcmc for one hour ended

			double pct_acc = 100.0 * n_acc/(n_acc+n_rej);
			printf("%3.0f%% ", pct_acc);

			//masud: the following block is used to adjust prop_std but Jinduan kept it commented. Going to uncomment it.
			/*
			if (iter < 5) {
				if (pct_acc < 20) prop_std /= 1.5;
				if (pct_acc > 60) prop_std *= 1.5;
			}
			*/
			
			
			

			t_pct_acc += pct_acc;
			//Pop_calc(pop);
			//Pop_report(pop);
			Pop_writeout(pop, estep_outfile, NULL, mcmc_chain_size, t);

			//update demand estimates with sample means
			Pop_Err perr = Pop_mean(pop, &buf[(vm->lb+t)*dim], dim, burn_in);
			/*
			for (int i=0; i<dim; ++i) {
				printf("%3.1f ", buf[(vm->lb+t)*dim + i]);
			}
			printf("\n");
			*/
			if (perr) return 301;
			
		} // End of loop: for (int t = 0; t < ws; ++t)
		
		//masud: adding the following line to adjust proposal standard deviation within EM-Iteration if percent acceptance are not within reasonable limits
		t_pct_acc /= ws;
		if (iter > 3) {
			if (t_pct_acc < 20.0) prop_std /= 1.5;
			if (t_pct_acc > 50.0) prop_std *= 1.5;
		}
		fprintf(fp162, "EM-Iteration: %d <--> prop_std = %f\n", iter, prop_std);
		
		//SM.-----------END OF E-Step---------------------------------------------------///////
		//SM. End of E step... MCMC population generation


		//M-step: re-estimate demand model
		//VarModel_Err vmerr = VarModel_estimate(vm, &buf[vm->lb*dim], dim, ws);
		
		double diff;
		double mdll; //demand-likelihood
		VarModel_Err vmerr = VarModel_estimate(vm, buf, dim, buf_size, &diff, &mdll);
		//memcpy(xd0, &buf[vm->lb*dim], sizeof(double)*dim*ws);
		if (vmerr) return 401;
		
		printf("\nVAR Model updated in M step (Error code: %d). Iteration %d\n",vmerr, iter);
		printf("\n");
		//VarModel_dump(vm);
		
		//masud: for log demand
		if (is_log_demand) {
			printf("exp(mu):\n");
			for (int i=0; i<dim; ++i) {
				printf("%6.3f ", exp(vm->mu[i]));
			}
		}
		double mhll = 0; // hydraulic likelihood
		
						 
		//sm./ tag:bug -- the program tries to use the var multiplier but it hasn't been initialized.
		// the following block produces the numbers in the console... if the format %d, %f, %f...
		for (int t = 0; t < ws; ++t) 
		{
			double mhllt = 0;
			 
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
					//int ind_sm = (vm->lb + t)*dim + i;
					//printf("t=%d, real_demand[%d] = %f,buf[%d]=%f\n",t,i,real_demand[i],ind_sm,buf[ind_sm]);
					real_demand[i] = buf[(vm->lb + t)*dim + i];
				}
			}
	
			//sm./ sv solves the network hydraulics with the 
			sv->run(real_demand, dim, &scada[t*n_ch], n_ch);
			sv->logL(&scada[t*n_ch], n_ch, &mhllt, true);	//masud: logL produces the %d %f %f output to console. 
			mhll += mhllt;
		}

		printf("Sum of changes in parameters: "
			"%.2f, dmll: %.2f, mhll: %.2f, total l:%.2f\n", 
			diff, mdll, mhll, mdll+mhll); 
		//if (t_pct_acc > 60) {
			//prop_std *= 1.3;
		//} else if (t_pct_acc < 15) {
		//	prop_std /= 1.5;
		//}
		//if (iter >= 5) {
			//prop_std /= 1.4; // reduce std dev of proposal for MCMC
		//}
		
	}
	// end EM cycle
	//sm. End of : for (int iter = 0; iter < n_iter; ++iter)
	fclose(fp162);

	print_mat(buf, dim, buf_size, "Buffer");

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
	delete sv;
	delete ds;
	// TODO: Network destructor
	//delete net;

	iniparser_freedict(dict);

	return 0;
}


