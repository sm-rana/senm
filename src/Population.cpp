#include <stdlib.h>
#include <tchar.h>
#include <math.h>
#include "Population.h"

Pop_Err Pop_new(int M_in, 
				int d1_in, 
				Population** pop_out) {
	return Pop_new(M_in, d1_in, 1, pop_out, 0);
}

Pop_Err Pop_new(int M_in, 
				int d1_in, 
				int d2_in, 
				Population** pop_out) {
	return Pop_new(M_in, d1_in, d2_in, pop_out, 0);
}

Pop_Err Pop_new(int M, 
				int d1, 
				int d2, 
				Population** pop_out, 
				int shadow) {
	Pop_Err err;

	Population* prec = (Population*)calloc(1, sizeof(Population));
	if (prec==NULL) {err = POP_MEM_NOT_ALLOCED; goto END;}

	if (!M || !d1 || !d2 || !pop_out) {err = POP_INPUT_ERR; goto END;}
	prec->M = M;
	prec->d1 = d1;
	prec->d2 = d2;


	int iw = 0;	// masud: treat iw as 0
	for (; iw<N_WORKERS; ++iw) {
		if (!shadow) prec->h[iw] = (double*) calloc(M*d1*d2, sizeof(double)); //masud: h[0][1000*7]
		prec->mean[iw] = (double*) calloc(d1*d2, sizeof(double)); //sm. basically mean[0][7]
		prec->sdev[iw] = (double*) calloc(d1*d2, sizeof(double)); //sm. basically sdev[0][7]
		if ((!shadow && !prec->h[iw]) || !prec->mean[iw] || !prec->sdev[iw]) // catching exceptions..ignore
		{err = POP_MEM_NOT_ALLOCED; goto END;}
		prec->p[iw] = 0;
		prec->_ps[iw] = 0;
	}
	
	prec->tmean = (double*) calloc(d1*d2, sizeof(double));	// tmean[7]
	prec->tsdev = (double*) calloc(d1*d2, sizeof(double));  // tsdev[7]
	if ( !prec->tmean || !prec->tsdev) 
	{err = POP_MEM_NOT_ALLOCED; goto END;}

	/*
	if (d2 == 1) {// vector series
		for (iw=0; iw<N_WORKERS; ++iw) {
			for (int iac=0; iac<POP_MAX_ACF; ++iac) {
				prec->acor[iw][iac] = (double*) calloc(d1, sizeof(double)); 
				if (!prec->acor[iw][iac]) {err = POP_MEM_NOT_ALLOCED; goto END;}
			}
		}
	}
	*/

	for (iw=0; iw<N_WORKERS; ++iw) {
		prec->clU1[iw] = (double*) calloc(d1*d2, sizeof(double));
		prec->clL1[iw] = (double*) calloc(d1*d2, sizeof(double));
		prec->clU2[iw] = (double*) calloc(d1*d2, sizeof(double));
		prec->clL2[iw] = (double*) calloc(d1*d2, sizeof(double));
		if (!prec->clU1[iw] || !prec->clL1[iw] || 
			!prec->clU2[iw] || !prec->clL2[iw])
		{err = POP_MEM_NOT_ALLOCED; goto END;}
	}
	prec->tclU1 = (double*) calloc(d1*d2, sizeof(double));
	prec->tclL1 = (double*) calloc(d1*d2, sizeof(double));
	prec->tclU2 = (double*) calloc(d1*d2, sizeof(double));
	prec->tclL2 = (double*) calloc(d1*d2, sizeof(double));
	if (!prec->tclU1 || !prec->tclL1 || 
		!prec->tclU2 || !prec->tclL2){
			err = POP_MEM_NOT_ALLOCED; goto END;}



	// init order stats
	if (!shadow) for (iw=0; iw<N_WORKERS; ++iw) {
		prec->_os[iw] = (int*)calloc(M*d1*d2, sizeof(int)); //os[][1000*7]
		if (!prec->_os[iw]) {
			err = POP_MEM_NOT_ALLOCED; goto END;}
		for (int im=0; im<M; ++im)
			for (int id=0; id<d1*d2; ++id)
				POPOS(prec, iw, im, id) = im; //pre-loaded indices for chains
				//masud: POPOS(prec, iw, im, id) = im EQUAL TO _os[0][0]=[0][1]=0; _os[0][2]=[0][3]=1; so on...
	}

	if (!shadow) {
		prec->_tos = (int*)calloc(M*d1*d2*N_WORKERS, sizeof(int));
		if (!prec->_tos) {
			err = POP_MEM_NOT_ALLOCED; goto END;
		}
		InitializeSRWLock(&prec->statLoc);
	}


	*pop_out = prec;
	return POP_OK;

END:
	if (prec!=NULL) free(prec);
	*pop_out = NULL;
	return err;

}

// chain order statistics ranking comparison, only used by Pop_chainCmp()
// and Pop_TchainCmp()
static int Pop_curChain;
static int Pop_curDim;
static Population* Pop_curInstance;

//---- internal functions ----
// comparison within a chain
static int Pop_chainCmp(const void* p1, const void* p2);

// comparison in total pop
static int Pop_TchainCmp(const void* p1, const void* p2);

static int partition(Population* pop, int w, int id, int l, int r) {
   int i, j, t;
   int pivot_im = POPOS(pop, w, l, id); 
   double pivot = POPH(pop, w, pivot_im, id);
   i = l; j = r+1;
		
   while( 1)
   {
   	do ++i; while( i <= r && POPH(pop, w, POPOS(pop, w, i, id), id) <= pivot);
   	do --j; while( POPH(pop, w, POPOS(pop, w, j, id), id) > pivot );
   	if( i >= j ) break;
   	t = POPOS(pop, w, i, id); 
	POPOS(pop, w, i, id) = POPOS(pop, w, j, id); 
	POPOS(pop, w, j, id) = t;
   }
   	t = POPOS(pop, w, l, id); 
	POPOS(pop, w, l, id) = POPOS(pop, w, j, id); 
	POPOS(pop, w, j, id) = t;
 
   return j;
}
static void quick_sort(Population* pop, int w, int id, int l, int r)
{
   int j;

   if( l < r ) 
   {
   	// divide and conquer
        j = partition(pop, w, id, l, r);
       quick_sort(pop, w, id, l, j-1);
       quick_sort(pop, w, id, j+1, r);
   }
	
}




void Pop_calc(Population* pop) {
	if (pop==NULL) return;

	int ps0[N_WORKERS]; //previous reporting pointer
	int ps1[N_WORKERS]; //current reporting pointer

	Population* sdpop;
	if (Pop_new(pop->M, pop->d1, pop->d2, &sdpop, 1) != POP_OK)  return;

	pop->perU1 = 1-(1-POP_CL_1)/2;  // percentiles
	pop->perU2 = 1-(1-POP_CL_2)/2;
	pop->perL2 = (1-POP_CL_2)/2;
	pop->perL1 = (1-POP_CL_1)/2;

	//calc single chain statistics
	for (int iw=0; iw<N_WORKERS; ++iw){
		ps0[iw] = pop->_ps[iw]; //previous reporting counter
		ps1[iw] = pop->p[iw]; //current reporting counter

		for (int id=0; id<pop->d1*pop->d2; ++id) { //iterate elements
			double mean0 = pop->mean[iw][id];
			double sum = mean0*ps0[iw];
			double sumvar =  POP_SQR(pop->sdev[iw][id])*ps0[iw]; //sum of variance
			double sumac[POP_MAX_ACF]; //sum of auto correlations

			/*
			if (pop->d2 == 1) for (int iac = 0; iac<POP_MAX_ACF; ++iac)
				sumac[iac] = (pop->acor[iw][iac][id])* ps0[iw];
				*/

			for (int im=ps0[iw]; im<ps1[iw]; ++im) { //iterate chain items
				sum += POPH(pop, iw, im, id);
				sumvar += POP_SQR( POPH(pop, iw, im, id)-mean0);

				if (pop->d2==1) for (int iac=0; iac<POP_MAX_ACF; ++iac) {
					if (im>iac) sumac[iac]+= 
						(POPH(pop,iw,im,id)-mean0)*(POPH(pop,iw,im-iac,id)-mean0);
				}
			}

			double mean1=sum/ps1[iw];

			//quick sort single chain os
			// TODO: make this part thread safe
			//Pop_curChain = iw;
			//Pop_curDim = id;
			//Pop_curInstance = pop;
			if (ps1[iw]>1) quick_sort(pop, iw, id, 0, ps1[iw]-1) ;

			// update single-chain order statistics to shadow pop
			int osL1 = int(ps1[iw]*pop->perL1);
			int osL2 = int(ps1[iw]*pop->perL2);
			int osU1 = int(ps1[iw]*pop->perU1);
			int osU2 = int(ps1[iw]*pop->perU2);
			sdpop->clU1[iw][id] = POPH(pop,iw,POPOS(pop,iw,osU1,id),id);
			sdpop->clU2[iw][id] = POPH(pop,iw,POPOS(pop,iw,osU2,id),id);
			sdpop->clL2[iw][id] = POPH(pop,iw,POPOS(pop,iw,osL2,id),id);
			sdpop->clL1[iw][id] = POPH(pop,iw,POPOS(pop,iw,osL1,id),id);

			//update mean, stddev, and auto-correlations to shadow pop
			sdpop->mean[iw][id] = mean1;
			sdpop->sdev[iw][id] = sqrt((sumvar - ps1[iw]*POP_SQR(mean0-mean1))/ps1[iw]);
			/*
			if (pop->d2 == 1) for (int iac=0; iac<POP_MAX_ACF; ++iac) {
                if (ps1[iw] > iac)
				sdpop->acor[iw][iac][id] =  
					(sumac[iac]-POP_SQR(mean0-mean1)*(ps1[iw]))/(ps1[iw]) ;
			}
			*/

		}

	}

	//TODO: chain total stats
	//calc total statistics
	for (int id = 0; id<pop->d1*pop->d2; ++id) {
		double tsum = 0;
		int pssum = 0;
		int tlen = N_WORKERS * pop->M; // length of total pop order stats

		int itos = 0; //iterator of total population order stats
		for (int iw=0; iw<N_WORKERS; ++iw)  {
			tsum += sdpop->mean[iw][id]*ps1[iw];
			pssum += ps1[iw];
			for (int im=ps0[iw]; im<ps1[iw]; ++im)  
				pop->_tos[id * tlen + (itos++)] = 
				im * N_WORKERS + iw; //coding in the position inside a chain and the chain no.
		}
		double tmean1 = tsum/pssum;

		double tsumvar = 0;
		for (int iw=0; iw<N_WORKERS; ++iw) {
			tsumvar += ( POP_SQR(sdpop->sdev[iw][id])*ps1[iw] + 
				ps1[iw]*POP_SQR(tmean1-sdpop->mean[iw][id]));
		}

		sdpop->tmean[id] = tmean1;
		sdpop->tsdev[id]= sqrt(tsumvar/pssum);

		//update order stat of total pop into shadow pop
		/*
		Pop_curInstance = pop;
		Pop_curDim = id;
		qsort ((void*) (& pop->_tos[id*tlen]), itos, sizeof(int), Pop_TchainCmp);

        int code;
        code = pop->_tos[id*tlen + int(itos*perU1)];
		sdpop->tclU1[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + int(itos*perU2)];
		sdpop->tclU2[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + int(itos*perL2)];
		sdpop->tclL2[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + int(itos*perL1)];
		sdpop->tclL1[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);
		*/
	}

	//update real pop stat using shadow pop - protected by a rwlock
	
	AcquireSRWLockExclusive(&pop->statLoc);
	for (int iw=0; iw<N_WORKERS; ++iw) {
		pop->_ps[iw] = ps1[iw]; // update report pointer
	}
	
	for (int id=0; id<pop->d1*pop->d2; ++id) {
		pop->tmean[id] = sdpop->tmean[id];
		pop->tsdev[id] = sdpop->tsdev[id];
		pop->tclU1[id] = sdpop->tclU1[id];
		pop->tclU2[id] = sdpop->tclU2[id];
		pop->tclL2[id] = sdpop->tclL2[id];
		pop->tclL1[id] = sdpop->tclL1[id];
		for (int iw=0; iw<N_WORKERS; ++iw) {
			pop->mean[iw][id] = sdpop->mean[iw][id];
			pop->sdev[iw][id] = sdpop->sdev[iw][id];
			pop->clU1[iw][id] = sdpop->clU1[iw][id];
			pop->clU2[iw][id] = sdpop->clU2[iw][id];
			pop->clL2[iw][id] = sdpop->clL2[iw][id];
			pop->clL1[iw][id] = sdpop->clL1[iw][id];
			/*
			if (pop->d2==1) for (int iac=0; iac<POP_MAX_ACF; ++iac) {
				pop->acor[iw][iac][id] = sdpop->acor[iw][iac][id];
			}
			*/
		}
	}
	ReleaseSRWLockExclusive(&pop->statLoc);

}

int Pop_TchainCmp(const void* p1, const void* p2) {
	int code1 = *(int*)p1;
	int code2 = *(int*)p2;
	int chain1 =  code1 % N_WORKERS;
	int pos1 = code1 / N_WORKERS;
	int chain2 = code2 % N_WORKERS;
	int pos2 = code2 / N_WORKERS;

	double x1 = POPH(Pop_curInstance, chain1, pos1, Pop_curDim);
	double x2 = POPH(Pop_curInstance, chain2, pos2, Pop_curDim);

	return (x1<x2?-1:(x1>x2?1:0));
}

/*
int Pop_chainCmp(const void* p1, const void* p2) {
	double x1 = POPH(Pop_curInstance, Pop_curChain, *(int*)p1, Pop_curDim);
	double x2 = POPH(Pop_curInstance, Pop_curChain, *(int*)p2, Pop_curDim);

	return (x1<x2?-1:(x1>x2?1:0));
}
*/

void Pop_report(Population* pop, FILE* target) {
	_ftprintf(target, 
		TEXT("\nPopulation %p Summary: \n")
		TEXT("Number of Chains/workers: [%u], Max chain size: %u, Dim: (%u*%u)\n"),
		pop, N_WORKERS, pop->M, pop->d1, pop->d2);

	int iw=0;
	int id=0;
	AcquireSRWLockShared(&pop->statLoc);
	for (iw=0; iw< N_WORKERS; ++iw) {
		_ftprintf(target, 
			TEXT("Chain[%u]: Counter %u/%u/%u (Reported/Computed/Max)\n"),
			iw, pop->_ps[iw], pop->p[iw], pop->M);
	}

    //mean, stddev, credible limits for total and single chains
    for (id=0; id < pop->d1*pop->d2; ++id) {
        _ftprintf(target,
			TEXT("\nDim (%u,%u): Mean-+Std.dev. {%6.2f-+%6.2f}"),
            id/pop->d1+1, id%pop->d1+1, pop->tmean[id], pop->tsdev[id]);

        //for (iw=0; iw<N_WORKERS; ++iw)
        //    _ftprintf(target, TEXT(", %6.2f-+%6.2f"), pop->mean[iw][id], pop->sdev[iw][id]);

        //_ftprintf(target,
        //    TEXT("\n\tP%4.2f={%6.2f,%6.2f}"), POP_CL_1, pop->tclU1[id], pop->tclL1[id]);

        for (iw=0; iw<N_WORKERS; ++iw)
            _ftprintf(target, TEXT(", (%d-%d-%d-%d%%)=[%6.2f, %6.2f, %6.2f, %6.2f]"),
			int(pop->perL1*100), int(pop->perL2*100), 
			int(pop->perU2*100), int(pop->perU1*100),
			pop->clL1[iw][id], pop->clL2[iw][id], pop->clU2[iw][id], pop->clU1[iw][id]);

        //_ftprintf(target,
        //    TEXT("\n\tP%4.2f={%6.2f,%6.2f}"), POP_CL_2, pop->tclU2[id], pop->tclL2[id]);

        //for (iw=0; iw<N_WORKERS; ++iw)
         //   _ftprintf(target, TEXT(", [%6.2f,%6.2f]"), );
	}
    _ftprintf(target, TEXT("\n"));

    // auto correlations
    int iac = 0;
    if (pop->d2 == 1) for (id=0; id< pop->d1; ++id) {
        _ftprintf(target, TEXT("Dim(%u) ACF:"), id);
        for (iac = 0; iac < POP_REP_ACF; ++iac) {
			_ftprintf(target, TEXT(" <%u>"), iac);
			/*
			for (iw=0; iw< N_WORKERS; ++iw) 
				_ftprintf(target, TEXT("%6.2f,"), 
				  pop->acor[iw][iac][id]/POP_SQR(pop->sdev[iw][id]));
				  */
		}
        _ftprintf(target, TEXT("\n"));
	}

	ReleaseSRWLockShared(&pop->statLoc);
}

Pop_Err Pop_mean(Population* pop, double* mean_out, int dim, int burn_in) {
	if (dim != pop->d1*pop->d2) return POP_DIM_MISMATCH;
	if (burn_in >= pop->M) return POP_CHAIN_TOO_SHORT;

	for (int iw=0; iw<N_WORKERS; ++iw) {
		if (pop->p[iw] != pop->M) return POP_MEAN_NOT_READY;
	}
	for (int id=0; id<dim; ++id) {
		mean_out[id] = 0;
		for (int im = burn_in; im<pop->M; ++im) {
			for (int iw = 0; iw < N_WORKERS; ++iw) {
				mean_out[id] += POPH(pop, iw, im, id);
			}
		}
		mean_out[id] /= N_WORKERS * (pop->M - burn_in);
	}
		
	return POP_OK;
}


void Pop_reset(Population* pop) {
	//reset all report pointers and stats
	for (int iw = 0; iw < N_WORKERS; ++iw) {
		pop->p[iw] = 0;
		pop->_ps[iw] = 0;
		memset(pop->mean[iw], 0, sizeof(double)*pop->d1*pop->d2); //masud: mean[7] is set to zero
		memset(pop->sdev[iw], 0, sizeof(double)*pop->d1*pop->d2); //masud: sdev[7] is set to zero
		memset(pop->h[iw], 0, sizeof(double)*pop->M*pop->d1*pop->d2); //masud: h[7000] is set to zero
		for (int im = 0; im < pop->M; ++im) {
			for (int id = 0; id < pop->d1*pop->d2; ++id) {
				POPOS(pop, iw, im, id) = im; //#define POPOS(pop, iw, im, id)   (pop->_os[iw][(id) + (im)*pop->d1*pop->d2])
				//masud: basically _os[0:6] =1; _os[7:13] = 2;...
			}
		}
	}

	memset(pop->tmean, 0, sizeof(double)*pop->d1*pop->d2);
	memset(pop->tsdev, 0, sizeof(double)*pop->d1*pop->d2);
	memset(pop->tclL1, 0, sizeof(double)*pop->d1*pop->d2);
	memset(pop->tclL2, 0, sizeof(double)*pop->d1*pop->d2);
	memset(pop->tclU1, 0, sizeof(double)*pop->d1*pop->d2);
	memset(pop->tclU2, 0, sizeof(double)*pop->d1*pop->d2);

	//TODO: all chains order stat



}

//masud: tag:outputfile
void Pop_writeout(Population* pop, 
				  char* filename, char* os_outfilename,
				  int chain_head, int timestep) {

	char out_full_filename[MAX_FILE_NAME_SIZE];
	sprintf(out_full_filename, "%s%d", filename, timestep);
	FILE* outfile = fopen(out_full_filename, "w");
	if (outfile == NULL) { return; }

	for (int im = 0; im < chain_head; ++ im) {
		for (int id = 0; id < pop->d1 * pop->d2; ++id) {
			fprintf(outfile, "%8.3f, ", POPH(pop, 0, im, id));
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);

	FILE* os_outfile; 
	char os_full_filename[MAX_FILE_NAME_SIZE];
	sprintf(os_full_filename, "%s%d", os_outfilename, timestep);
	if (os_outfilename != NULL)  {
		os_outfile = fopen(os_full_filename, "w");
		for (int id = 0; id < pop->d1 * pop->d2; ++id) {
			for (int im = 0; im < chain_head; ++im) {
				fprintf(os_outfile, "%8.3f ", 
					POPH(pop, 0, POPOS(pop, 0, im, id), id));
			}
			fprintf(os_outfile, "\n");
		}
		fclose(os_outfile);
	}
}

void Pop_writeout2(Population* pop,
	char* filename, char* os_outfilename,
	int chain_head, int timestep, int iter) {

	char out_full_filename[MAX_FILE_NAME_SIZE];
	sprintf(out_full_filename, "%s_%d_%d",filename,iter, timestep);
	FILE* outfile = fopen(out_full_filename, "w");
	if (outfile == NULL) { return; }

	for (int im = 0; im < chain_head; ++im) {
		for (int id = 0; id < pop->d1 * pop->d2; ++id) {
			fprintf(outfile, "%8.3f ", POPH(pop, 0, im, id));
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
	
}


void Pop_del(Population** ppop) {
    Population* pop = (*ppop);
	int iw=0;
	for (; iw<N_WORKERS; ++iw) {
		free(pop->h[iw]);
		free(pop->mean[iw]);
		free(pop->sdev[iw]);

		/*
		if (pop->d2==1) {
			for (int iac=0; iac<POP_MAX_ACF; ++iac)
				free(pop->acor[iw][iac]);
		} 
		*/

		free(pop->clU1[iw]);
		free(pop->clU2[iw]);
		free(pop->clL1[iw]);
		free(pop->clL2[iw]);

		free(pop->_os[iw]);

	}

	free(pop->tmean);
	free(pop->tsdev);
	free(pop->_tos);

	free(pop->tclL1);
	free(pop->tclU1);
	free(pop->tclL2);
	free(pop->tclU2);

    (*ppop) = NULL;

}
