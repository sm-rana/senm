#include <stdlib.h>
#include <tchar.h>
#include <math.h>
#include "Population.h"

Pop_Err Pop_new(unsigned M_in, unsigned d1_in, Population** pop_out) {
	return Pop_new(M_in, d1_in, 1, pop_out, 0);
}
Pop_Err Pop_new(unsigned M_in, unsigned d1_in, unsigned d2_in, Population** pop_out) {
	return Pop_new(M_in, d1_in, d2_in, pop_out, 0);
}
Pop_Err Pop_new(unsigned M, unsigned d1, unsigned d2, Population** pop_out, int shadow) {
	Pop_Err err;

	Population* prec = (Population*)calloc(1, sizeof(Population));
	if (prec==NULL) {err = POP_MEM_NOT_ALLOCED; goto END;}

	if (!M || !d1 || !d2 || !pop_out) {err = POP_INPUT_ERR; goto END;}
	prec->M = M;
	prec->d1 = d1;
	prec->d2 = d2;


	unsigned iw = 0;
	for (; iw<N_WORKERS; ++iw) {
		if (!shadow) prec->h[iw] = (double*) calloc(M*d1*d2, sizeof(double));
		prec->mean[iw] = (double*) calloc(d1*d2, sizeof(double));
		prec->sdev[iw] = (double*) calloc(d1*d2, sizeof(double));
		if ((!shadow && !prec->h[iw]) || !prec->mean[iw] || !prec->sdev[iw])
		{err = POP_MEM_NOT_ALLOCED; goto END;}
		prec->p[iw] = 0;
		prec->_ps[iw] = 0;
	}

	prec->tmean = (double*) calloc(d1*d2, sizeof(double));
	prec->tsdev = (double*) calloc(d1*d2, sizeof(double));
	if ( !prec->tmean || !prec->tsdev) 
	{err = POP_MEM_NOT_ALLOCED; goto END;}

	if (d2 == 1) {// vector series
		for (iw=0; iw<N_WORKERS; ++iw) {
			for (unsigned iac=0; iac<POP_MAX_ACF; ++iac) {
				prec->acor[iw][iac] = (double*) calloc(d1, sizeof(double)); 
				if (!prec->acor[iw][iac]) {err = POP_MEM_NOT_ALLOCED; goto END;}
			}
		}
	}

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
		prec->_os[iw] = (unsigned*)calloc(M*d1*d2, sizeof(unsigned));
		if (!prec->_os[iw]) {
			err = POP_MEM_NOT_ALLOCED; goto END;}
		for (unsigned im=0; im<M; ++im)
			for (unsigned id=0; id<d1*d2; ++id)
				POPOS(prec, iw, im, id) = im; //pre-loaded indices for chains
	}

	if (!shadow) {
		prec->_tos = (unsigned*)calloc(M*d1*d2*N_WORKERS, sizeof(unsigned));
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
static unsigned Pop_curChain;
static unsigned Pop_curDim;
static Population* Pop_curInstance;

void Pop_calc(Population* pop) {
	if (pop==NULL) return;

	unsigned ps0[N_WORKERS]; //previous reporting pointer
	unsigned ps1[N_WORKERS]; //current reporting pointer

	Population* sdpop;
	if (Pop_new(pop->M, pop->d1, pop->d2, &sdpop, 1) != POP_OK)  return;

	double perU1 = 1-(1-POP_CL_1)/2;  // percentiles
	double perU2 = 1-(1-POP_CL_2)/2;
	double perL2 = (1-POP_CL_2)/2;
	double perL1 = (1-POP_CL_1)/2;

	//calc single chain statistics
	for (unsigned iw=0; iw<N_WORKERS; ++iw){
		ps0[iw] = pop->_ps[iw]; //previous reporting counter
		ps1[iw] = pop->p[iw]; //current reporting counter

		for (unsigned id=0; id<pop->d1*pop->d2; ++id) { //iterate elements
			double mean0 = pop->mean[iw][id];
			double sum = mean0*ps0[iw];
			double sumvar =  POP_SQR(pop->sdev[iw][id])*ps0[iw]; //sum of variance
			double sumac[POP_MAX_ACF]; //sum of auto correlations

			if (pop->d2 == 1) for (unsigned iac = 0; iac<POP_MAX_ACF; ++iac)
				sumac[iac] = (pop->acor[iw][iac][id])* ps0[iw];

			for (unsigned im=ps0[iw]; im<ps1[iw]; ++im) { //iterate chain items
				sum += POPH(pop, iw, im, id);
				sumvar += POP_SQR( POPH(pop, iw, im, id)-mean0);

				if (pop->d2==1) for (unsigned iac=0; iac<POP_MAX_ACF; ++iac) {
					if (im>iac) sumac[iac]+= 
						(POPH(pop,iw,im,id)-mean0)*(POPH(pop,iw,im-iac,id)-mean0);
				}
			}

			double mean1=sum/ps1[iw];

			//quick sort single chain os
			Pop_curChain = iw;
			Pop_curDim = id;
			Pop_curInstance = pop;
			qsort( (void*)(& POPOS(pop,iw,0,id)), ps1[iw], sizeof(unsigned), Pop_chainCmp);  

			// update single-chain order statistics to shadow pop
			sdpop->clU1[iw][id] = POPH(pop,iw,POPOS(pop,iw,unsigned(ps1[iw]*perU1),id),id);
			sdpop->clU2[iw][id] = POPH(pop,iw,POPOS(pop,iw,unsigned(ps1[iw]*perU2),id),id);
			sdpop->clL2[iw][id] = POPH(pop,iw,POPOS(pop,iw,unsigned(ps1[iw]*perL2),id),id);
			sdpop->clL1[iw][id] = POPH(pop,iw,POPOS(pop,iw,unsigned(ps1[iw]*perL1),id),id);

			//update mean, stddev, and auto-correlations to shadow pop
			sdpop->mean[iw][id] = mean1;
			sdpop->sdev[iw][id] = sqrt((sumvar - ps1[iw]*POP_SQR(mean0-mean1))/ps1[iw]);
			if (pop->d2 == 1) for (unsigned iac=0; iac<POP_MAX_ACF; ++iac) {
                if (ps1[iw] > iac)
				sdpop->acor[iw][iac][id] =  
					(sumac[iac]-POP_SQR(mean0-mean1)*(ps1[iw]))/(ps1[iw]) ;
			}

		}

	}

	//calc total statistics

	for (unsigned id = 0; id<pop->d1*pop->d2; ++id) {
		double tsum = 0;
		unsigned pssum = 0;
		unsigned tlen = N_WORKERS * pop->M; // length of total pop order stats

		unsigned itos = 0; //iterator of total population order stats
		for (unsigned iw=0; iw<N_WORKERS; ++iw)  {
			tsum += sdpop->mean[iw][id]*ps1[iw];
			pssum += ps1[iw];
			for (unsigned im=ps0[iw]; im<ps1[iw]; ++im)  
				pop->_tos[id * tlen + (itos++)] = 
				im * N_WORKERS + iw; //coding in the position inside a chain and the chain no.
		}
		double tmean1 = tsum/pssum;

		double tsumvar = 0;
		for (unsigned iw=0; iw<N_WORKERS; ++iw) {
			tsumvar += ( POP_SQR(sdpop->sdev[iw][id])*ps1[iw] + 
				ps1[iw]*POP_SQR(tmean1-sdpop->mean[iw][id]));
		}

		sdpop->tmean[id] = tmean1;
		sdpop->tsdev[id]= sqrt(tsumvar/pssum);

		//update order stat of total pop into shadow pop

		Pop_curInstance = pop;
		Pop_curDim = id;
		qsort ((void*) (& pop->_tos[id*tlen]), itos, sizeof(unsigned), Pop_TchainCmp);

        unsigned code;
        code = pop->_tos[id*tlen + unsigned(itos*perU1)];
		sdpop->tclU1[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + unsigned(itos*perU2)];
		sdpop->tclU2[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + unsigned(itos*perL2)];
		sdpop->tclL2[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);

        code = pop->_tos[id*tlen + unsigned(itos*perL1)];
		sdpop->tclL1[id] = POPH(pop, code%N_WORKERS, code/N_WORKERS, id);
	}

	//update real pop stat using shadow pop - protected by a rwlock
	AcquireSRWLockExclusive(&pop->statLoc);
	for (unsigned iw=0; iw<N_WORKERS; ++iw) {
		pop->_ps[iw] = ps1[iw]; // report pointer
	}
	for (unsigned id=0; id<pop->d1*pop->d2; ++id) {
		pop->tmean[id] = sdpop->tmean[id];
		pop->tsdev[id] = sdpop->tsdev[id];
		pop->tclU1[id] = sdpop->tclU1[id];
		pop->tclU2[id] = sdpop->tclU2[id];
		pop->tclL2[id] = sdpop->tclL2[id];
		pop->tclL1[id] = sdpop->tclL1[id];
		for (unsigned iw=0; iw<N_WORKERS; ++iw) {
			pop->mean[iw][id] = sdpop->mean[iw][id];
			pop->sdev[iw][id] = sdpop->sdev[iw][id];
			pop->clU1[iw][id] = sdpop->clU1[iw][id];
			pop->clU2[iw][id] = sdpop->clU2[iw][id];
			pop->clL2[iw][id] = sdpop->clL2[iw][id];
			pop->clL1[iw][id] = sdpop->clL1[iw][id];
			if (pop->d2==1) for (unsigned iac=0; iac<POP_MAX_ACF; ++iac) {
				pop->acor[iw][iac][id] = sdpop->acor[iw][iac][id];
			}
		}
	}
	ReleaseSRWLockExclusive(&pop->statLoc);
}

int Pop_TchainCmp(const void* p1, const void* p2) {
	unsigned code1 = *(unsigned*)p1;
	unsigned code2 = *(unsigned*)p2;
	unsigned chain1 =  code1 % N_WORKERS;
	unsigned pos1 = code1 / N_WORKERS;
	unsigned chain2 = code2 % N_WORKERS;
	unsigned pos2 = code2 / N_WORKERS;

	double x1 = POPH(Pop_curInstance, chain1, pos1, Pop_curDim);
	double x2 = POPH(Pop_curInstance, chain2, pos2, Pop_curDim);

	return (x1<x2?-1:(x1>x2?1:0));
}


int Pop_chainCmp(const void* p1, const void* p2) {
	double x1 = POPH(Pop_curInstance, Pop_curChain, *(unsigned*)p1, Pop_curDim);
	double x2 = POPH(Pop_curInstance, Pop_curChain, *(unsigned*)p2, Pop_curDim);

	return (x1<x2?-1:(x1>x2?1:0));
}

void Pop_report(Population* pop, FILE* target) {
	_ftprintf(target, 
		TEXT("\nPopulation %p Summary: \n")
		TEXT("Number of Chains/workers: [%u], Max chain size: %u, Dim: (%u*%u)\n"),
		pop, N_WORKERS, pop->M, pop->d1, pop->d2);

	unsigned iw=0;
	unsigned id=0;
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

        for (iw=0; iw<N_WORKERS; ++iw)
            _ftprintf(target, TEXT(", %6.2f-+%6.2f"), pop->mean[iw][id], pop->sdev[iw][id]);

        _ftprintf(target,
            TEXT("\n\tP%4.2f={%6.2f,%6.2f}"), POP_CL_1, pop->tclU1[id], pop->tclL1[id]);

        for (iw=0; iw<N_WORKERS; ++iw)
            _ftprintf(target, TEXT(", [%6.2f,%6.2f]"), pop->clU1[iw][id], pop->clL1[iw][id]);

        _ftprintf(target,
            TEXT("\n\tP%4.2f={%6.2f,%6.2f}"), POP_CL_2, pop->tclU2[id], pop->tclL2[id]);

        for (iw=0; iw<N_WORKERS; ++iw)
            _ftprintf(target, TEXT(", [%6.2f,%6.2f]"), pop->clU2[iw][id], pop->clL2[iw][id]);
	}
    _ftprintf(target, TEXT("\n"));

    // auto correlations
    unsigned iac = 0;
    if (pop->d2 == 1) for (id=0; id< pop->d1; ++id) {
        _ftprintf(target, TEXT("Dim(%u) ACF:"), id);
        for (iac = 0; iac < POP_REP_ACF; ++iac) {
			_ftprintf(target, TEXT(" <%u>"), iac);
			for (iw=0; iw< N_WORKERS; ++iw) 
				_ftprintf(target, TEXT("%6.2f,"), pop->acor[iw][iac][id]);
		}
        _ftprintf(target, TEXT("\n"));
	}

	ReleaseSRWLockShared(&pop->statLoc);
}

void Pop_del(Population* pop) {
	unsigned iw=0;
	for (; iw<N_WORKERS; ++iw) {
		free(pop->h[iw]);
		free(pop->mean[iw]);
		free(pop->sdev[iw]);

		if (pop->d2==1) {
			for (unsigned iac=0; iac<POP_MAX_ACF; ++iac)
				free(pop->acor[iw][iac]);
		} 

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

}
