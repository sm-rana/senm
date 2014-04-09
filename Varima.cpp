#include <cmath>
#include <cfloat>
#include "Varima.h"
#include "clapack.h"



Varima::Err Varima::setPara(
	double* phi_in, // n*p 
	double* theta_in, // n*q
	double* cov_in, // n*n
	double* mean_in // n*1
	)
{

	if (s==1 && !cov_in ) 
		// some parameters are not set
			return PARA_NOT_SET; 
	if (!(n>0)) 
		return NUM_VARIABLE_ERR;

	/* allocate memory */
	if (p>0) phi = (double*)calloc( n * (p+1), sizeof(double));
	if (q>0) theta = (double*)calloc( n * (q+1) , sizeof(double));
	if (d>0) dwts=(int*) calloc( s*d+1, sizeof(int));

	mean = (double*)calloc( n, sizeof(double));

	if (s==1) {
		cov = (double*)calloc(n*n, sizeof(double));
		lcov = (double*)calloc(n*n, sizeof(double));
	}

	ncma = s*q+f+1; //number of columnes for ma
	ncmw = s*p+f+1;
	ncmz = s*d+f+1;
	ma = (double*)calloc(n*(s*q+f+1), sizeof(double));
	if (s==1) tpa = (double*)calloc(n, sizeof(double)); //inner-most only
	mw = (double*)calloc(n*(s*p+f+1), sizeof(double));
	mz = (double*)calloc(n*(s*d+f+1), sizeof(double));

	if (!(mean && ma && mw && mz)) return CANT_ALLOC_MEM;
	if (p>0 && !phi) return CANT_ALLOC_MEM;
	if (q>0 && !theta) return CANT_ALLOC_MEM;
	if (d>0 && !dwts) return CANT_ALLOC_MEM;
	if (s==1 && (!cov || !lcov)) return CANT_ALLOC_MEM;

	int i=0;
	/* copy parameters */
	for (i=0; i<n; ++i) { //iterate through each single time series
		if (p>0) {
			phi[i * (p+1)] = 1;  //phi_0
			for (int j=1; j<=p; ++j) {
				phi[i * (p+1)  +  j] = -phi_in[i * p + j-1]; //
			}
		}
		if (q>0)  {
			theta[ i * (q+1)] = 1;  //theta_0
			for (int j=1; j<=q; ++j) {
				theta[i * (q+1) + j] = -theta_in[i * q + j-1];
			}
		}
	}
	if (mean_in) for (i=0; i<n; ++i) mean[i]=mean_in[i];
	int j;

	maxCov = minCov = 0;
	if (s==1) {
		maxCov = -DBL_MAX; minCov = DBL_MAX;
		for (i=0; i<n; ++i) //row iter
			for (j=0; j<=i; ++j) {//col iter
				double tp = cov_in[i*n +j];
				cov[i*n + j] = tp; //row-major
				lcov[j*n + i] = tp; //col-major
				if (tp > maxCov) maxCov = tp ;
				if (tp < minCov) minCov = tp;
			}
	}

	/* set dwts */
	if (d>0) for (i=0; i<=d; ++i) dwts[i*s]=(int)nchoosek(d, i)*(i%2 ? -1 : 1);

	/* Cholesky decomposition of cov, gotoblas/clapack is 
	column-major */
	integer info;
	if (s==1) {
		dpotrf_("L", (integer*)&n, lcov, (integer*)&n, &info);

		if (info) return COV_NOT_POSITIVE_DEF;
	}



	return OK;

}


Varima::Err Varima::setPara(double* phi_in, double* theta_in, double* cov_in) {
	return setPara(phi_in, theta_in, cov_in, NULL);
}

Varima::~Varima() {
	free(phi); 
	free(theta);
	free(cov);
	free(lcov);
	free(mean);
	free(ma);
	free(tpa);
	free(mw);
	free(mz);
}

Varima::Err Varima::initStateG(Rvgs* rg_in, double** z0, int size) {
	// z0 is a pointer to rows
	if (rg_in == NULL) return INVALID_PTR;
	// set random number generator for the inner most model
	Varima* inner = this;
	while (inner->s!=1) inner=inner->in;
	inner->rg = rg_in;

	//if ( z0 == NULL) return OK; // don't init

	int j; //dimension iterator
	int k; //time iterator
	int kc;  //column iterator

	if (z0) // has inital matrix
		for (j=0; j<n; ++j)  // fill initial value with very last values
			for (k=ncmz-1, kc=size-1; k>=0 && kc>=0; --k, --kc) 
				mz[j*ncmz + k] = z0[j][kc];


	for (m=0; m < INIT_CYCLES * (p+1) * (q+1) * (d+1);) 
		if (generate(NULL)) return GEN_ERROR;

	return OK;
}

//Varima::Err Varima::initStateLogG(Rvgs* rg_in, double** z0, int size) {
//	// z0 is a pointer to rows
//	if (rg_in == NULL) return INVALID_PTR;

//	// set random number generator for the inner most model
//	Varima* inner = this;
//	while (inner->s!=1) inner=inner->in;
//	inner->rg = rg_in;

//	//if ( z0 == NULL) return OK; // don't init

//	int j; //dimension iterator
//	int k; //time iterator
//	int kc;  //column iterator

//	if (z0) // has inital matrix
//		for (j=0; j<n; ++j)  // fill initial value with very last values
//			for (k=ncmz-1, kc=size-1; k>=0 && kc>=0; --k, --kc) {
//				if (z0[j][kc] <= 0) 
//					mz[j*ncmz + k] = LOG_OF_ZERO;
//				else
//					mz[j*ncmz + k] = log(z0[j][kc]);
//			};


//    for (m=0; m < INIT_CYCLES * (p+1) * (q+1) * (d+1);) 
//        if (generate(NULL)) return GEN_ERROR;

//    return OK;
//}


Varima::Err  Varima::generate(double* z_out) {

	int j; //dimension iterator
	int k; //parameter iterator

	//advance cyclic buffer pointers
	pma++; pma %= ncma;
	pmw++; pmw %= ncmw;
	pmz++; pmz %= ncmz;
	++m; //incresee counter

	// update ma
	if (s==1) { //inner-most model, no seasonality
		//generate current a
		for (j=0; j<n; ++j) tpa[j]=rg->Normal(0,1); //get iid ran.vec
		doublereal alpha(1), beta(0); 
		integer incx=1;
		// use cholesky lower-triangular matrix to generate correlated vector white noise
		dgemv_("N", (integer*)&n, (integer*)&n, &alpha, lcov, (integer*)&n,
			tpa, &incx, &beta, ma, (integer*)&ncma); // L matrix * vector
	} else {// use embeded model
		if (in->generate(NULL)) return GEN_ERROR;
		for (j=0; j<n; ++j) 
			ma[j*ncma + pma]=(in->mz)[j*in->ncmz + in->pmz]; 
	}


	// update mw
	for (j=0; j<n; ++j) {
		double tpw = ma[j*ncma + pma]; // current a
		for (k=1; k<=p; ++k) { //AR parameters
			tpw -= 
				mw[j*ncmw + (pmw - s*k + ncmw)%ncmw]  //use cyclic mem get previous w
			* phi[j*(p+1) + k];
		}
		for (k=1; k<=q; ++k) { //MA parameters
			tpw += ma[j*ncma + (pma - s*k + ncma)%ncma]
			* theta[j*(q+1) + k];
		}
		mw[j*ncmw + pmw] = tpw;
	}

	//update mz and fill return array
	for (j=0; j<n; ++j) {
		double tpz = mw[j*ncmw + pmw];
		for (k=1; k<=d; ++k)  // differencing operator
			tpz -= mz[j*ncmz + (pmz - s*k + ncmz)%ncmz] * dwts[k*s];

		mz[j*ncmz + pmz] = tpz;
		if (z_out != NULL) z_out[j] = tpz + mean[j];
	}

	return OK;
}

Varima::Err Varima::generateExp(double* z_out) {
	Err err;
	if (err = generate(NULL)) return err;

	if (z_out == NULL) return INVALID_PTR;
	for (int j=0; j<n; ++j)

		z_out[j] = exp(mz[j*ncmz + pmz]); 

	return OK;
}


unsigned Varima::nchoosek(unsigned n, unsigned k) {
	// does not work for very large n and k

	unsigned long long res = 1;
	unsigned long long tp = 1;

	if (n == 0 || n < k) return 0;
	if (k > n/2 + 1) k = n-k;
	if (k == 0) return 1;

	while (k) {
		res *= (n--) ;
		tp *= (k--);
	}
	return (unsigned)(res/tp);
}
