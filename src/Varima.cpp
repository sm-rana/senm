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
	dwts=(int*) calloc( s*d+1, sizeof(int));
    dwts[0] = 1;

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
	if (!dwts) return CANT_ALLOC_MEM;
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

void Varima::reportEWI(Err err) {
	TCHAR ewiTxts[DUMMY_LAST][MAX_ERR_STRING_SIZE]= {
        TEXT("OK."),
        TEXT("Parameters have not been set."),
        TEXT("Dimension of the model must be a positive integer."),
        TEXT("Could not allocate memory."),
        TEXT("The covariance matrix is not positively definitive, try using two independent models."),
        TEXT("The pointer is invalid."),
        TEXT("Could not generate a vector."),
        TEXT("Not enough data to initialize the model."),
        TEXT("Log of zero when initializating."),
	};

	TCHAR txt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];
    _stprintf(txt, TEXT("<ARI>Varima:%s"), ewiTxts[err]);

    ewi(txt, threadN);
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

Varima::Err Varima::initStateG(Rvgs* rg_in, double** z0, int size, int bLog) {
    if (z0 == NULL) { //set random number generator and exit
        Varima* inner = this;
        while (inner->s!=1) inner=inner->in;
        inner->rg = rg_in;
		return OK;
	}
    double* tempZ0 = (double*)calloc(n*size, sizeof(double));
    for (int i =0; i<n; ++i) 
        for (int j = 0; j < size; ++j) {
            if (bLog) { 
                if (z0[i][j] <= 0) return LOG_INIT_WITH_ZERO;
				else tempZ0[i*size + j] = log(z0[i][j]);
			} else {
                tempZ0[i*size + j] = z0[i][j];
			}
		}
    Err err = initStateG(rg_in, tempZ0, size);
    free(tempZ0);
    return err;
}

Varima::Err Varima::initStateG(Rvgs* rg_in, double* z0, int size) {
    // This initialization method introduces errors for MA models
	// z0 is a pointer to rows
	if (rg_in == NULL) return INVALID_PTR;

	// set random number generator for the inner most model
    if (in == NULL) rg = rg_in; //inner-most model

	//if ( z0 == NULL) return OK; // don't init

    if (size < s*d + s*p + 1) return INIT_SIZE_TOO_SMALL;


	int j; //dimension iterator
	int k; //time iterator - varima buffer
	int kc;  //column iterator - source data
    int h; // linear combination iterator
    // temporary storage of intermedite init data
    int nctpmw = size-s*d;
    int nctpma = size-s*d-s*p;
    double* tpmw = (double*)calloc(n*nctpmw, sizeof(double));
    double* tpma = (double*)calloc(n*nctpma, sizeof(double));

    for (j = 0; j<n; ++j) {
        // fill z-buffer
        for (k = ncmz-1, kc=size-1; k>= 0; --k, --kc) 
            mz[j*ncmz+k] = z0[j*size+kc];

        // fill w-buffer
        for (k=1; k<= size-s*d; ++k) {
            double tpw = 0;
            for (h=0; h<=s*d; ++h) tpw += dwts[h]*z0[j*size+size-k-h];
            tpmw[j*nctpmw + nctpmw-k] = tpw;
            if (ncmw-k>=0) mw[j*ncmw+ncmw-k] = tpw;
		}

        // fill a-buffer
        for (k=1; k<= size-s*d-s*p; ++k) {
            double tpa = 0;
            for (h=0; h<=p; ++h) tpa += phi[h]*tpmw[j*nctpmw+nctpmw-k-h*s];
            tpma[j*nctpma + nctpma-k] = tpa;
            if (ncma-k>=0) ma[j*ncma+ncma-k] = tpa;
		}
	}

    Err err;
	if (in!=NULL) err = in->initStateG(rg_in, tpma, nctpma);
    free(tpmw);
    free(tpma);

//	if (z0) // has inital matrix
//		for (j=0; j<n; ++j)  // fill initial value with very last values
//			for (k=ncmz-1, kc=size-1; k>=0 && kc>=0; --k, --kc) {
//                if (bLog == 0)
//                    mz[j*ncmz + k] = z0[j][kc];
//				else if (z0[j][kc] <= 0)  // log(z0) init
//					mz[j*ncmz + k] = LOG_OF_ZERO;
//				else
//					mz[j*ncmz + k] = log(z0[j][kc]);
//			}

//	for (m=0; m < INIT_CYCLES * (p+1) * (q+1) * (d+1);) 
//		if (generate(NULL)) return GEN_ERROR;

	return OK;
}


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
		// use lcov matrix to generate correlated vector white noise
		dgemv_("N", (integer*)&n, (integer*)&n, &alpha, lcov, (integer*)&n,
			tpa, &incx, &beta, &ma[pma], (integer*)&ncma); // L matrix * vector
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
