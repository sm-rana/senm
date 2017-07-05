/* SenM - Sensor-driven hydraulic Model of water distribution systems
Copyright (C) 2013  Jinduan Chen <jinduan.uc@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License v2
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html)
for more details.   */


#pragma once
#include <Windows.h>
#include "stdio.h"
#include "SenmCoreIncs.h"

//#define POPH2(pop, iw, im, row, col)   (pop->h[iw][((row)*pop->d1+(col))*pop->M + (im)])
//#define POPH(pop, iw, im, id)   (pop->h[iw][(id)*pop->M + (im)])
//#define POPOS(pop, iw, im, id)   (pop->_os[iw][(id)*pop->M + (im)])


/// two macros for easier matrix/vector indexing
#define POPH2(pop, iw, im, row, col)   (pop->h[iw][((col)*pop->d1+(row)) + (im)*pop->d1*pop->d2])
#define POPH(pop, iw, im, id)   (pop->h[iw][(id) + (im)*pop->d1*pop->d2])
#define POPOS(pop, iw, im, id)   (pop->_os[iw][(id) + (im)*pop->d1*pop->d2])


// default max-lag of autocovariance computation (hours in a month)
//#define POP_MAX_ACF (24*7*4)  
#define POP_MAX_ACF 5  
// report acf
#define POP_REP_ACF 5

// Credible limits probability range
#define POP_CL_1  0.95
#define POP_CL_2  0.50

//SQR and CUB
#define POP_SQR(x) ((x)*(x))
#define POP_CUB(x) ((x)*(x)*(x))

/// Data storage and statistic reporting of vector/matrix chains 
struct Population {

	int M;  ///> max length of the chain;
	int d1;  ///> size of 1st dimension
	int d2;  ///> size of 2nd dimension, =1 if vector chain

	double * h[N_WORKERS]; ///> head of vector/matrix chains, each with Dim d1*d2*M
	int p[N_WORKERS]; ///> chain pointer, how many elements in h has been produced

    /// report pointer, corresponds to the stats below in each chain
	int _ps[N_WORKERS]; 

    /// statistics access lock, reading following stats must acquire the shared lock first
	SRWLOCK statLoc;

	// statistics
	double* mean[N_WORKERS]; ///> chain means, dim = d1*d2
	double* tmean; ///> total means (of all chains)
	double* sdev[N_WORKERS];  ///> chain std deviation, dim = d1*d2
	double* tsdev; ///> total std deviation (of all chains)

    /// autocovariances for vector time series, only available 
	/// for 1d data (d2==1)
	//double* acor[N_WORKERS][POP_MAX_ACF]; 

    /// order statistics - ranks in a single chain 
    /** _os[iw][id  + im*d1*d2] is the index (im in h[]) of the im-th smallest 
	 value in h[] for the id-th element in the iw-th chain */
	int *_os[N_WORKERS];  

    /// order statistics - ranks in all chains
    int *_tos;

	double perU1, perU2, perL1, perL2;

	double* clU1[N_WORKERS];  ///> credible limit upper: 1st. each chain (i.e., 97.5% percentile)
	double* clL1[N_WORKERS];  ///> credible limit lower: 1st
	double* clU2[N_WORKERS];  ///> credible limit upper: 2nd
	double* clL2[N_WORKERS];  ///> credible limit lower: 2nd

	double* tclL1;  ///> credible limit lower: 1st. total 
	double* tclU2;  ///> credible limit upper: 2nd. total 
	double* tclL2;  ///> credible limit lower: 2nd. total 
	double* tclU1;  ///> credible limit upper 1st. total

};

///> Error code for struct Population
enum Pop_Err {
		POP_OK,
		POP_MEM_NOT_ALLOCED,
		POP_INPUT_ERR,
		POP_DIM_MISMATCH,
		POP_CHAIN_TOO_SHORT,
		POP_MEAN_NOT_READY,

		POP_DUMMY_LAST};


/// create Population - Factory method
/** M_in, d1_in, d2_in must >0 */
Pop_Err Pop_new(int M_in, int d1_in, int d2_in, Population** pop_out);


/// create Population - Factory method, 
/** if shadow = 1, create a shadow instance (no h[], os[], tos data storage, 
   only stats storage, used by Pop_calc() */
Pop_Err Pop_new(int M_in, int d1_in, int d2_in, Population** pop_out, int shadow);
Pop_Err Pop_new(int M_in, int d1_in, Population** pop_out) ;

void Pop_del(Population** pop);

//report means excluding samples for burning-in
Pop_Err Pop_mean(Population* pop, double* mean_out, int dim, int burn_in=0);

///> report statistics
void Pop_report(Population* pop, FILE* target=stdout);

//TODO: update also stats pointers p[]
///> Update all statistics. 
void Pop_calc(Population* pop);

///> Reset all statistics
void Pop_reset(Population* pop);


///> Dump the poplulation to a data file
void Pop_writeout(Population* pop, 
				  char* filename, char* os_outfilename=NULL,
				  int interval=0, int t=0);

