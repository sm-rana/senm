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
	
	int is_log_demand = iniparser_getint(dict, "XDemands:is_log_demand", 0);
	int is_dev_demand = iniparser_getint(dict, "XDemands:is_dev_demand", 0);


	int n_iter = iniparser_getint(dict, "EMMC:em_iterations", 0);
	if (n_iter <= 0) {
		return 106;
	}


	int mcmc_chain_size = iniparser_getint(dict, "EMMC:max_chain_size", 0);
	if (mcmc_chain_size <= 0) {
		return 101;
	}

	int burn_in = iniparser_getint(dict, "EMMC:burn_in", 0);
	if (burn_in <= 0 || burn_in > mcmc_chain_size) {
		return 107;
	}

/*
	int init_xd_size = iniparser_getint(dict, 
		"EMMC:initial_demand_data_length", 0);
	if (init_xd_size <= 0) {
		return 102;
	}
	*/

	int est_win_size = iniparser_getint(dict, 
		"EMMC:estimation_win_size", 0);
	if (est_win_size <= 0) {
		return 104;
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

	if (mode == 0) net->report();

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

	int pattern_size;
	double *multiplier;
	double* real_demand; // real demands for a single time step
	if (is_dev_demand) {
		fscanf(xmf, "%d", &pattern_size);
		multiplier = (double*)calloc(pattern_size, sizeof(double));
		real_demand = (double*)calloc(dim, sizeof(double));
		for (int i=0; i<pattern_size; ++i) {
			fscanf(xmf, "%lf", &multiplier[i]);
		}
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

	if (mode == 0) {
		printf("SenM checked out successfully.");
		goto END;
	}

	// Loading scada panel data

	printf("[0]Loading scada panel data from the channels for estimation...\n");
	printf("  Time interval %d s, Current time (Year month day) in EM: %s. Using "
		   "previous %d snapshots.\n", 
		ds->dt,  est_end_time, est_win_size);
	Tstamp tcur;  // estimation end time for EM, usually the current time stamp
	sscanf(est_end_time, "%d %d %d %d %d %d", 
		&tcur.year, &tcur.month, &tcur.day, 
		&tcur.hour, &tcur.minute, &tcur.second);
	tcur.fraction = 0;
	
	double* scada;
	int n_ch = ds->n_chan;
	int ws = est_win_size;
	scada = (double*)calloc(n_ch * ws, sizeof(double));
	dsec = ds->fillSnapshots(tcur, ws, scada);


	//demand data buffer for var model look-back and estimation
	int buf_size = vm->lb + ws;
	double* buf = (double*)calloc(buf_size*dim, sizeof(double));

	//load init demand - for a continuously running program this may
	// be the previous demand estimates
	memcpy(buf, &xd0[(est_win_size - vm->lb)*dim], 
				sizeof(double)*vm->lb*dim);
	memcpy(&buf[vm->lb*dim], xd0, sizeof(double)*ws*dim);


	//init random number generator
	Rvgs* rng = new Rvgs(2);

	// EM cycle
	for (int iter = 0; iter < n_iter; ++iter) {
		printf("EM iteration %d: \n", iter);
		double t_pct_acc = 0;

		// E-step: sampling demands
		for (int t = 0; t < ws; ++t) {

			//first sample
			Pop_reset(pop);
			//memcpy(&POPH(pop, 0, 0, 0), mu , sizeof(double)*dim);
			memcpy(&POPH(pop, 0, 0, 0), &buf[(vm->lb+t)*dim], 
				sizeof(double)*dim);

			// MCMC
			double cur_ll = 0;  //current log-likelihood
			double last_ll = -DBL_MAX;  //ll from the previous sample
			int n_acc = 1;
			int n_rej = 0;
			for (int im = 0; im < mcmc_chain_size ; ++im) {

				double* cur_xd = &POPH(pop, 0, im, 0);
				int has_wrong_demand = 0;

				if (is_dev_demand) {
					for (int i=0; i<dim; ++i) {
						real_demand[i] = multiplier[t%pattern_size]*mu[i] + cur_xd[i];
					}
				}

				if (is_dev_demand) {
					for (int id = 0; id < dim; ++id) { //reject problematic demands
						double tp;
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

				if (! has_wrong_demand) {
					// hydraulic likelihood
					if (is_dev_demand) {
						sv->run(real_demand, dim, &scada[t*n_ch], n_ch);
					} else if (is_log_demand) {
						sv->runlogd(cur_xd, dim, &scada[t*n_ch], n_ch);
					} else {
						sv->run(cur_xd, dim, &scada[t*n_ch], n_ch);
					}
					sv->logL(&scada[t*n_ch], n_ch, &cur_ll);

					//compute demand likelihood
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
				if (im == mcmc_chain_size - 1) break;
				else {//propose a new demand panel
					for (int id=0; id<dim; ++id) {
						POPH(pop, 0, im+1, id) = POPH(pop, 0, im, id) 
							+ rng->Normal(0, prop_std);
					}
				}

			}
			double pct_acc = 100.0 * n_acc/(n_acc+n_rej);
			printf("%3.0f%% ", pct_acc);

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
			
		}
		t_pct_acc /= ws;

		//M-step: re-estimate demand model
		//VarModel_Err vmerr = VarModel_estimate(vm, &buf[vm->lb*dim], dim, ws);
		double diff;
		double mdll; //demand-likelihood
		VarModel_Err vmerr = VarModel_estimate(vm, buf, dim, buf_size, &diff, &mdll);
		//memcpy(xd0, &buf[vm->lb*dim], sizeof(double)*dim*ws);
		if (vmerr) return 401;
		//printf("VAR Model updated in M step. Iteration %d\n", iter);
		printf("\n");
		//VarModel_dump(vm);
		if (is_log_demand) {
			printf("exp(mu):\n");
			for (int i=0; i<dim; ++i) {
				printf("%6.3f ", exp(vm->mu[i]));
			}
		}
		double mhll = 0; // hydraulic likelihood

		for (int t = 0; t < ws; ++t) {
			double mhllt = 0;
			for (int i=0; i<dim; ++i) {
				real_demand[i] = multiplier[t%pattern_size]*mu[i] 
					+ buf[(vm->lb+t)*dim+i];
			}
	
			sv->run(real_demand, dim, &scada[t*n_ch], n_ch);
			sv->logL(&scada[t*n_ch], n_ch, &mhllt, true);
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
			prop_std /= 1.4; // reduce std dev of proposal for MCMC
		//}
	}

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


