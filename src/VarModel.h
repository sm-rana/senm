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
struct VarModel {
	int n; ///< dimensions
	int ns; ///< number of seasonal levels
	int *s; ///< length of seasonal periods
	int *p; ///< number of autoregressive parameters on each level

	double* phi;  ///< array of all AR matrices on all seasonal levels
	///< phi[0] = I
	double* cov;  ///< covariance matrix (column-major)
	double* mu;  /// mean of the VAR

	// derived variables

	int phi_cnt; ///< cnt of phi array
	int lb; ///< Look-back length
	double* lcov; ///< covariance matrix - Cholesky decomposition (bottom-left)
	double sqrt_d; ///< sqrt of cov determinant
	double* poly;  ///< characteristic polynomial. only non-zeros are stored 
	int* pos; ///< look-back position of ar parameters
	int npos; ///<number of elements in pos and poly
	///< the pos and para arrays provide a more concise representation of 
	///< the poly array
};

///Macros for matrix indexing
///#define SARI_H(sari, i, j) ((sari)->H[(i)+(j)*3])

enum VarModel_Err {
	VARMODEL_OK,
	VARMODEL_ALLOC_ERROR,
	VARMODEL_PHI_ARRAY_ERROR,
	VARMODEL_COV_MATRIX_ERROR,
	VARMODEL_ATRIX_ERROR,
	VARMODEL_ARRAY_ERROR,
	VARMODEL_NOT_ENOUGH_SAMPLE,
	VARMODEL_DIM_ERR,
	VARMODEL_NO_AR_PARA
};

/// Create a Sari model
VarModel_Err VarModel_new(int n_in, int ns_in, 
				  int* s_in, int* p_in, 
				  VarModel** varm_out);

/// Destroy a Sari model
void VarModel_del(VarModel** varm);

/// Set Autogressive parameters
VarModel_Err VarModel_setPhi(VarModel* varm, double* phi_in, int phi_size);

/// Set covariance matrix and compute cholesky decomposition
VarModel_Err VarModel_setCov(VarModel* varm, double* cov_in, int cov_len);

/// Estimate the VAR mean to the average of panel data/or the param if dim=1
VarModel_Err VarModel_estMu(VarModel* varm, double* panel, int panel_dim,
							int panel_len);

/// Compute log-likelihood, thread-safe
VarModel_Err VarModel_logLikelihood(VarModel* varm, double* panel, 
				int panel_dim, int panel_len, 
				double* logl_out, double* a_out=NULL, int n_a_out=0 /*effective no. w. noise*/);

/// Estimate parameters from the provided data. original parameters poly and cov
/// are overwritten
// diff - Euclidean distance between old and new parameters
VarModel_Err VarModel_estimate(VarModel* varm, double* panel, 
				int panel_dim, int panel_size, double* diff, double* mdll);

/*
/// polynomial multiplication
void VarModel_poly_mlpl(int* poly_a, ///< a polynomial (phi indices)
					  int* len_a,  ///< length of poly_a
					  int* poly_b, ///< another polynomial (phi idx)
									///< the first "1" is omitted
					  int len_b, ///<length of poly_b
					  int s ///< step of poly_b
					  ) ;
					  */

/// Dump the data of a VarModel struct
void VarModel_dump(VarModel* varm);

//Sari_Err VarModel_ErrRpt(Sari* sari, 

void print_mat(double* p_mat, ///<matrix to be printed, column-major
					  int n_row, 
					  int n_col,
					  char* name);
	
