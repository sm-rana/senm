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

#include <stdlib.h>
#include <stdio.h>
#include <rvgs.h>
#include "SenMCoreIncs.h"

/// this number determines the number of initialization points
/// = INIT_CYCLES * (p+1) * (q+1) * (d+1)
#define INIT_CYCLES  5 
#define LOG_OF_ZERO  -1000000  /// so exp(LOG_OF_ZERO) ~ 0

/// Multi-dimensional seasonal auto-regressive integrated moving average model
/** The model has the following usage: 
    1. generate time series with given coefficients 
        (i.e., temporal/spatial correlations)
        (e.g., used by a SCADAgenerator)
    2. forecast the series with given parameters
      ( as in the inference step of SenM)

   "look-back" memory of previous observations and white noise realizations
   are stored for AR models

    3. Determine the likelihood of a demand realization
      ( as in the E(MCMC)-step of SenM, and M-step MLE computation) 
    3 uses provided state info. 		*/
class Varima {
protected:
	int n;  ///> dimensionality
	int p;  ///> number of autoregressive parameters (AR)
	int q;  ///> number of moving-average parameters (MA)

    /** seasonality. if s > 1, this is an outer seasonal model 
	  and there should be an inner model pointed by in */
	int s;  
	Varima* in;  ///> inner model, models can be nested recursively
	int d;  ///> degree of integration (I)

	int f;  ///> Max forecasting horizon

	/** matrix of AR parameters, dim=(n*(p+1)), row-major,
        AR parameters are stored as left-hand-side values (negative of the input) */
	double* phi;  

	/** matrix of MA parameters, dim=(n*(q+1)), row-major,
	  MA parameters are stored as right-hand-side values (negative of the input) */
	double* theta; 

	int*	dwts;  ///> d-weights, expansion of the polynomial (1-B)^d, dim = 1*(d+1)

    /// covaraince matrix for the white noise, dim=(n*n), row-major
      /** since it is a harmitian positive definitive matrix,
       only the bottom-left diagonal half is used.
       applicable only for inner-most model */
	double* cov;  
	double maxCov, minCov; ///>max/min of var/cov
    
    
    /// cholesky decomposition of cov, dim=(n*n), row-major
    /** applicable only for inner-most model and generator
      only the bottom-left diagnal half is used. */
	double* lcov;  

    /// mean or initial "baseline" for the model, dim=(n*1)
    /** applicable only for generators, normally mean is only
       used in outmost model */
	double* mean;  

	/* time series state info*/
	double* ma;  ///>cyclic buffer/memory of a_t, dim=(n* (s*q+f+1) )
	double* tpa;  ///> temporary a_t storage, dim=(n * 1)
	double* mw;  ///>cyclic memory of w_t, dim=(n* (s*p+f+1) )
	double* mz;  ///>cyclic memory of z_t, dim=(n* (s*d+f+1) )
	int ncma, ncmw, ncmz;  ///>number of columes for ma, mw, mz
	
	int pma, pmw, pmz;  ///>pointers
	int m;  ///>data counter

    Rvgs* rg;  ///>random number generator

	enum Err {
		OK=0,
		PARA_NOT_SET,
		NUM_VARIABLE_ERR,
		CANT_ALLOC_MEM,
		COV_NOT_POSITIVE_DEF,
		INVALID_PTR,
		GEN_ERROR,
        INIT_SIZE_TOO_SMALL,
        LOG_INIT_WITH_ZERO,

        DUMMY_LAST
	};

    ///combinational number calculation/lookup
	unsigned nchoosek(unsigned n, unsigned k); 

private:
	Varima() {};  //can not instantiate w/o parameters
public:
    void reportEWI(Err err);
    unsigned threadN; ///>thread number, for ewi messages
	//constructors:  inner-most model
	Varima(int n_in, int p_in, int q_in):  //ARMA
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(0), f(0),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1), threadN(0),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {};

	Varima(int n_in, int p_in, int q_in, int d_in):  //ARIMA
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(d_in), f(0),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1), threadN(0),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {};

	Varima(int n_in, int p_in, int q_in, int d_in, int f_in):  //ARIMA forecaster
		n(n_in), p(p_in), q(q_in), in(NULL) /*inner most*/, d(d_in), f(f_in),
		m(0), pma(-1), pmw(-1), pmz(-1), 
		s(1), threadN(0),
	    phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL) {}; //

	///constructors: outer models, must use inner model's n and f
	Varima(Varima* inin, int p_in, int q_in, int s_in, int d_in): 
		n(inin->n), p(p_in), q(q_in), in(inin) /*inner most*/, d(d_in), f(inin->f),
		m(inin->m), pma(-1), pmw(-1), pmz(-1), 
		s(s_in), threadN(0),
		phi(NULL), theta(NULL), cov(NULL), lcov(NULL), mean(NULL),
		ma(NULL), mw(NULL), mz(NULL){};

	/// copy constructors
	Varima(const Varima& in, int f_in);  //generator->forecaster

    /// look-back window size
	int getLB() {return s*p+s*d+(in==NULL ?0:in->getLB());};

	/// get s, s>0 means a outer model
	int getS() {return s;}; 

	/// get f. f>0 means a forecaster
	int getF() {return f;};

    /// get dimensionality
	int getN() {return n;};

    ///get the total number of AR parameters including inner models
	int getNphi () { 
		if (in != NULL) 
			return p + in->getNphi();
		else 
			return p;
	};

    ///get the number of MA parameters
	int getNtheta() { 
		if (in != NULL) 
			return q + in->getNtheta();
		else 
			return q;
	}


	///get the i-th AR parameter of the k-th element,
	/** notice phi has p+1 columns, and i is valid only when >1
	   if i>p, then we return the inner model's (i-p)th parameters */
	double getPhi(int k, int i) {
		if (i<0 || k<0 || k>=n) return 0;
		if (i==0) return 1;
		if (i > p) {
			if (in == NULL) return 0;
			else return in->getPhi(k, i-p);
		} else {
			return -phi[k*(p+1) + i];
		}
	}

	/// get the i-th MA parameter of the k-th element
	double getTheta(int k, int i) {
		if (i<0 || k<0 || k>=n) return 0;
		if (i==0) return 1;
		if (i > q) {
			if (in == NULL) return 0;
			else return in->getTheta(k, i-q);
		} else {
			return -theta[k*(q+1) + i];
		}
	}

	///get the (i+1)-th linear parameter, including AR and MA parameter
	double getPara(int k, int i){
		if (i+1> this->getNphi())
			return this->getTheta(k, i+1-this->getNphi());
		else 
			return this->getPhi(k, i+1);
	}
    ///get level of the model, the inner-most model's level is 0
	int getLevel() {
		if (in == NULL) return 0;
		else return in->getLevel()+1;
	}
	///get string description of the (i+1)-th AR parameter, text must be pre-allocated
	void getPhiDesc(int i, char* desc) {
		int totalNphi = this->getNphi();
		if (i>= totalNphi) {
			sprintf(desc, "");
            return; //i-th parameter not exist
		}
		if (i<p) {
			sprintf(desc, "Level-%d, Phi_%d=", this->getLevel(), i+1);
		} else {
			if (in!=NULL) in->getPhiDesc(i-p, desc);
		}
	}
    ///get string description of the (i+1)-th MA parameter
	void getThetaDesc(int i, char* desc) {
		
		if (i>= this->getNtheta()) {
			sprintf(desc, "");
			return; //i-th parameter not exist
		}
		if (i<q) {
			sprintf(desc, "Level-%d, Theta_%d=", this->getLevel(), i+1);
		} else {
			if (in!=NULL) in->getThetaDesc(i-q, desc);
		}
	}
    ///get string description of the (i+1)-th parameter
	void getParaDesc(int i, char* desc) {
		int totalNphi = this->getNphi();
		if (i>=totalNphi) this->getThetaDesc(i-totalNphi, desc);
		else this->getPhiDesc(i, desc);
	}
    ///get the covariance(i,j)
	double getCov(int i, int j) {
		if (i>=0 && i<n && j>=0 && j<n) {
			if (in==NULL) 
				return cov[i*n + j];
			else
				return in->getCov(i,j);
		} else {
			return -1;
		}
	}

    /// @{
    /// Get the max/min value in the covariance matrix
	double getMaxCov() {return in?in->getMaxCov():maxCov;}
	double getMinCov() {return in?in->getMinCov():minCov;}
    /// @}

    /// @{
	/// set parameters, allocate memory
	Err setPara(double* phi_in, double* theta_in, double* cov_in, double* mean_in);
	Err setPara(double* phi_in, double* theta_in, double* cov_in);

    /// @}

	//  for both the generator and the forecaster

	/** initialize state info use z0, z0 must be a pointer array with n rows to 
	  'size' columns or NULL (init with 0's); */

	Err initStateLogG(Rvgs* rg_in, double** z0, int size) {
        return initStateG(rg_in, z0, size, 1);
	};  //log(z0) is used to init

	Err initStateG(Rvgs* rg_in, double** z0, int z0_col_size, int bLog=0);  

	Err initStateG(Rvgs* rg_in, double* z0, int z0_col_size);  

    /// time series data generation
	/** parameters of the models must be set already by setPara()
       and interal "look-back" memory initialized by initStateG()
    */
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
