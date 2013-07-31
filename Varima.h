#pragma once
//#include <stdlib.h>
#include <rvgs.h>

class Varima {
	/* Vector seasonal auto-regressive integrated moving average model*/
	/* The model has the following usage: 
		1. generate time series with given temporal/spatial correlations
		  ( as in a SCADAgenerator)
		2. forecast the series with given parameters
		  ( as in the inference step of SenM)
		1 and 2 uses self state memory of previous observations.

		3. Determine the likelihood of a demand realization
		  ( as in the E(MCMC)-step of SenM, and M-step MLE computation) 
		3 uses provided state info. 		*/

protected:
	int n;  // dimensionality
	int p;  // number of autoregressive parameters (AR)
	int q;  // number of moving-average parameters (MA)
	int s;  // seasonality. if s > 1, this is an outer seasonal model and 
			// there should be an inner model pointed by im
	Varima* in;  // inner model, models can be nested recursively
	int d;  // degree of integration (I)

	int f;  //forecasting length

	//AR parameters are stored as left-hand-side values (negative of the input)
	double* phi;  //matrix of AR parameters, dim=(n*(p+1)), row-major, 
	//MA parameters are stored as right-hand-side values (negative of the input)
	double* theta; //matrix of MA parameters, dim=(n*(q+1))

	int*	dwts;  // d-weights, expansion of (1-B)^d, dim = 1*(d+1)
	double* cov;  //covaraince matrix for the white noise, dim=(n*n)
				  // since it is a symetric positive definitive matrix,
				   //only the bottom-left diagonal half is used.
				  // applicable only for inner-most model
	double maxCov, minCov; //max/min of var/cov
	double* lcov;  // cholesky decomposition of cov, dim=(n*n)
					//applicable only for inner-most model and generator
					//only the bottom-left diagnal half is used.
	double* mean;  // mean or initial "baseline" for the model, dim=(n*1)
				   // applicable only for generators, mean is usually only
				  // used in outmost model

	/* time series state info*/
	double* ma;  //cyclic buffer/memory of a_t, dim=(n* (s*q+f+1) )
	double* tpa;  // temporary a_t storage, dim=(n * 1)
	double* mw;  //cyclic memory of w_t, dim=(n* (s*p+f+1) )
	double* mz;  //cyclic memory of z_t, dim=(n* (s*d+f+1) )
	int ncma, ncmw, ncmz;  //number of columes for ma, mw, mz
	
	int pma, pmw, pmz;  //pointers
	int m;  //counter

	// this number determines the number of initialization points
	// = INIT_CYCLES * (p+1) * (q+1) * (d+1)
	static const int INIT_CYCLES = 5; 

	Rvgs* rg;  //random number generator

	enum Err {
		OK=0,
		PARA_NOT_SET,
		NUM_VARIABLE_ERR,
		CANT_ALLOC_MEM,
		COV_NOT_POSITIVE_DEF,
		INVALID_PTR,
		GEN_ERROR
	
	};

	static const int LOG_OF_ZERO = -1000000;  //exp(LOG_OF_ZERO)~0

private:
	Varima() {};  //can not instantiate w/o parameters
	unsigned nchoosek(unsigned n, unsigned k); //combinational number
public:
	//constructors:  inner-most model
	Varima(int n_in, int p_in, int q_in):  //ARMA
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(0), f(0),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {};

	Varima(int n_in, int p_in, int q_in, int d_in):  //ARIMA
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(d_in), f(0),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {};

	Varima(int n_in, int p_in, int q_in, int d_in, int f_in):  //ARIMA forecaster
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(d_in), f(f_in),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {}; //

	//constructors: outer models, must use inner model's n and f
	Varima(Varima* inin, int p_in, int q_in, int s_in, int d_in): 
		n(inin->n), p(p_in), q(q_in), in(inin) /*inner most*/, d(d_in), f(inin->f),
		m(inin->m), pma(-1), pmw(-1), pmz(-1), 
		s(s_in),
		phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL){};

	// copy constructors
	Varima(const Varima& in, int f_in);  //generator->forecaster

	// get s, tell if it is a outer model
	int getS() {return s;}; 

	// get f. tell if it is a forecaster
	int getF() {return f;};

	int getN() {return n;};


	int getNphi () { //get the number of AR parameters including inner models
		if (in != NULL) 
			return p + in->getNphi();
		else 
			return p;
	};

	double getPhi(int k, int i) {
	//get the i-th AR parameter of the k-th element,
	// notice phi has p+1 columns, and i is valid only when >1
	// if i>p, then we return the inner model's (i-p)th parameters
		if (i<0 || k<0 || k>=n) return 0;
		if (i==0) return 1;
		if (i > p) {
			if (in == NULL) return 0;
			else return in->getPhi(k, i-p);
		} else {
			return -phi[k*(p+1) + i];
		}
	}

	int getNtheta() { //get the number of MA parameters
		if (in != NULL) 
			return q + in->getNtheta();
		else 
			return q;
	}

	double getTheta(int k, int i) {
	// get the i-th MA parameter of the k-th element
		if (i<0 || k<0 || k>=n) return 0;
		if (i==0) return 1;
		if (i > q) {
			if (in == NULL) return 0;
			else return in->getTheta(k, i-q);
		} else {
			return -theta[k*(q+1) + i];
		}
	}
	
	double getPara(int k, int i){//get the i-th parameter (i starts with 0)
		if (i+1> this->getNphi())
			return this->getTheta(k, i+1-this->getNphi());
		else 
			return this->getPhi(k, i+1);
	}

	int getLevel() {//get level of the model
		if (in == NULL) return 0;
		else return in->getLevel()+1;
	}
	
	void getPhiDesc(int i, char* desc) {//get string description of the i-th AR parameter
		int totalNphi = this->getNphi();
		if (i>= totalNphi) return; //i-th parameter not exist
		if (i<p) {
			sprintf(desc, "Level-%d, Phi_%d=", this->getLevel(), i+1);
		} else {
			if (in!=NULL) in->getPhiDesc(i-p, desc);
		}
	}

	void getThetaDesc(int i, char* desc) {//get string description of the i-th MA parameter
		
		if (i>= this->getNtheta()) return; //i-th parameter not exist
		if (i<q) {
			sprintf(desc, "Level-%d, Theta_%d=", this->getLevel(), i+1);
		} else {
			if (in!=NULL) in->getThetaDesc(i-q, desc);
		}
	}

	void getParaDesc(int i, char* desc) {//get string description of the i-th parameter
		int totalNphi = this->getNphi();
		if (i>=totalNphi) this->getThetaDesc(i-totalNphi, desc);
		else this->getPhiDesc(i, desc);
	}
		
	double getCov(int i, int j) {//get the covariance of i-j
		if (i>=0 && i<n && j>=0 && j<n) {
			if (in==NULL) 
				return cov[i*n + j];
			else
				return in->getCov(i,j);
		} else {
			return -1;
		}
	}

	double getMaxCov() {return in?in->getMaxCov():maxCov;}
	double getMinCov() {return in?in->getMinCov():minCov;}

	// set parameters
	Err setPara(double* phi_in, double* theta_in, double* cov_in, double* mean_in);
	Err setPara(double* phi_in, double* theta_in, double* cov_in);

	//  for generator and forecaster)
	// initialize state info use z0, z0 must be a matrix or pointer to pointer
	// with 'size' columns;
	//Err initStateG(Rvgs* rg_in, double* z0, int size);  
	Err initStateG(Rvgs* rg_in, double** z0, int size);  
	Err initStateLogG(Rvgs* rg_in, double** z0, int size);  //log(z0) is used to init


	// for generator
	Err generate(double* z_out);  // generate a vector 
	// and update state memory, z_out should have at least (n*1) space
	// if z_out is null, proceed anyway (internal mem is updated)

	Err generateExp(double* z_out);  //z_out = exp(z) for each element

	// for forecaster)
	double* getCB();  // get confidence bands. return array dim=(n*f) 
	double* getFC(double* obs);  // input new obs and get forecasts. 
	// input array dim=(n*1), return array dim=(n*f)

	// for likelihood)
	// get required raw data length of MLE
	int getLen();
	// get log-likelihood
	double logLikelihood(double*);



	~Varima();

};
