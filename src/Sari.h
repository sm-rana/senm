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


///Multi-dimensional multi-seasonal autoregressive models 
/** The same set of autogressive parameters (phi) are used for all dimensions.
	phi represents temporal correlations.
	A covariance matrix (cov) are used for spatial correlations. 
*/
struct Sari {
	unsigned n; ///< dimensions
	unsigned ns; ///< number of seasonal levels
	unsigned* H; ///< structural matrix (column-major), 3 lines are s, d, p
	double* phi;  ///< concatenated array of all AR parameters on all seasonal levels
	double* cov;  ///< covariance matrix (column-major)

	// derived variables

	unsigned phi_len; ///< length of phi array
	unsigned lb; ///< Look-back length
	double* lcov; ///< covariance matrix - Cholesky decomposition (bottom-left)
	double* poly;  ///< characteristic polynomial poly[0]=1
	unsigned* pos; ///< look-back position of ar parameters
	double* para; ///< the ar parameters at the positions specified in pos
	unsigned npos; ///<number of elements in pos and para
	///< the pos and para arrays provide a more concise representation of 
	///< the poly array
};

///Macros for matrix indexing
#define SARI_H(sari, i, j) ((sari)->H[(i)+(j)*3])

enum Sari_Err {
	SARI_OK,
	SARI_ALLOC_ERROR,
	SARI_H_MATRIX_ERROR,
	SARI_PHI_ARRAY_ERROR
};

/// Create a Sari model
Sari_Err Sari_new(unsigned n_in, unsigned ns_in, 
				  unsigned* H_in, unsigned H_len, Sari** model_out);

/// Destroy a Sari model
void Sari_del(Sari** sari);

/// Set Autogressive parameters
Sari_Err Sari_setPhi(Sari* sari, double* phi_in, unsigned phi_len);

/// Set covariance matrix and compute cholesky decomposition
Sari_Err Sari_setCov(Sari* sari, double* cov_in, unsigned cov_len);

/// Compute log-likelihood, thread-safe
Sari_Err Sari_logLikelihood(Sari* sari, double* sample, 
							unsigned sample_n, unsigned sample_len);


/// Dump the data of a Sari struct
void Sari_dump(Sari* sari);

//Sari_Err Sari_ErrRpt(Sari* sari, 


