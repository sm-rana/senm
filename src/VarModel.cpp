#include "limits.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "SenMCoreIncs.h"
#include "clapack_3.2.1\clapack.h"
#include "VarModel.h"


VarModel_Err VarModel_new(int n_in,  ///< dimension
				  int ns_in, ///< seasonality
				  int* s_in, ///
				  int* p_in,
				  VarModel** varm_out) {

	VarModel* p = (VarModel*)calloc(1, sizeof(VarModel));
	if (p == NULL) return VARMODEL_ALLOC_ERROR;

	p->n = n_in; p->ns = ns_in;
	p->s = (int*)calloc(ns_in + 1, sizeof(int));
	p->p = (int*)calloc(ns_in + 1, sizeof(int));
	if (p->s == NULL || p->p == NULL) 
		return VARMODEL_ALLOC_ERROR;

	// copy s, p
	memcpy(p->s, s_in, (ns_in + 1)*sizeof(int));
	memcpy(p->p, p_in, (ns_in + 1)*sizeof(int));

	// compute look-back length
	int lb = 0;
	int phi_cnt = 0;
	for (int i_ns = 0; i_ns <= p->ns; ++i_ns)  {
		// s * p 
		lb += p->s[i_ns] * p->p[i_ns];
		phi_cnt += p->p[i_ns];
	}
	p->lb = lb;
	p->phi_cnt = phi_cnt;

	// allocate mem for cov, phi, lcov
	p->cov = (double*)calloc(n_in * n_in, sizeof(double));
	p->lcov = (double*)calloc(n_in * n_in, sizeof(double));
	p->mu = (double*)calloc(n_in, sizeof(double));
	p->phi = (double*)calloc(n_in * n_in * phi_cnt, sizeof(double));
	p->poly = NULL; // don't know the size as now
	if (p->cov == NULL || p->lcov == NULL || p->phi == NULL || p->mu == NULL )
		return VARMODEL_ALLOC_ERROR;

	*varm_out = p;
	return VARMODEL_OK;
}

void VarModel_del(VarModel** varm) {
	free((*varm)->s);
	free((*varm)->p);

	free((*varm)->cov);
	free((*varm)->lcov);

	free((*varm)->poly);
	free((*varm)->phi);

	free(*varm);
}

VarModel_Err VarModel_setPhi(VarModel* varm, 
							 double* phi_in, 
							 int phi_size) {
	long n = (long)varm->n;
	int phi_cnt = varm->phi_cnt;
	if (phi_size != n * n * (phi_cnt)) 
		return VARMODEL_PHI_ARRAY_ERROR;
	memcpy(varm->phi, phi_in, sizeof(double)*phi_size);
	
	double* tp_poly; // temporary array of matrices
	int* tp_pos; // temporary array of non-zero para positions
	int n_items = 1;
	for (int i_ns = 0; i_ns <= varm->ns; ++ i_ns) {
		n_items *= (varm->p[i_ns] + 1);
	}

	tp_poly = (double*)calloc(n*n*n_items, sizeof(double));
	tp_pos = (int*)calloc(n_items, sizeof(int));
	varm->poly = (double*)calloc(n*n*n_items, sizeof(double));
	varm->pos = (int*)calloc(n_items, sizeof(int));
	for (int i = 0; i < n; ++i) {
		tp_poly[i + n*i] = 1;
	} // iden. mat 
	int len = 1; // length of the tp_pos (also no. matrices in tp_poly)

	double alpha = -1;
	double beta = 0;
	int ip0 = 0;
	for (int i_ns = 0; i_ns <= varm->ns; ++i_ns) {
		int p = varm->p[i_ns];
		int s = varm->s[i_ns];

		for (int ip =0; ip < p; ++ip) { // new polynomial
			for (int itp = 0; itp < len; ++itp) { 
				// matrix multiplication
				dgemm_("n", "n", &n, &n, &n, &alpha, 
					&tp_poly[itp*n*n], &n,  // left matrix
					&varm->phi[(ip + ip0)*n*n], &n,  // right matrix
					&beta, 
					&tp_poly[(itp+(ip+1)*len)*n*n], &n); //out
				tp_pos[itp+(ip+1)*len] = tp_pos[itp] + (ip+1)*s;
			}
		}
		ip0 += p; // update phi matrix pointer
		len *= (p+1); // roll flat tp_poly and tp_pos (vec() operator)
	}

	//bubble sort tp_pos, len should be small so little impact to performance
	int x;
	integer inc = 1;
	integer n2 = n*n;
	for (int i=len-1; i>=1; --i) {
		for (int j=0; j<i; ++j) {
			if (tp_pos[j] > tp_pos[j+1]) {
				x = tp_pos[j+1];
				tp_pos[j+1] = tp_pos[j];
				tp_pos[j] = x;
				dswap_(&n2, &tp_poly[n*n*j], &inc, 
							&tp_poly[n*n*(j+1)], &inc);
			}
		}
	}

	//merge items of same pos/order
	int p_poly = 0;
	int cur_order = 0;
	double a = 1;
	varm->pos[0] = 0;
	for (int i=0; i<len; ++i) {
		if (tp_pos[i] == cur_order) {
			// add to exising item
			daxpy_(&n2, &a, &tp_poly[i*n2], &inc, 
							&varm->poly[p_poly*n2], &inc);
		} else { // new item
			++ p_poly;
			cur_order = tp_pos[i];
			varm->pos[p_poly] = cur_order;
			memcpy(&varm->poly[p_poly*n2], &tp_poly[i*n2], 
				n2 * sizeof(double));
		}
	}

	varm->lb = cur_order;
	varm->npos = p_poly + 1;

	free(tp_poly);
	free(tp_pos);
	return VARMODEL_OK;
}

static void print_mat(int* p_mat, ///<matrix to be printed, column-major
					  int n_row, 
					  int n_col, 
					  char* name) {
	if (p_mat == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s matrix (positive integers): (%u x %u)\n", 
		name, n_row, n_col);
	for (int ir = 0; ir < n_row; ++ir) {
		for (int ic = 0; ic < n_col; ++ic) {
			printf("%4u ", p_mat[ir + n_row * ic]);
		}
		printf("\n");
	}
}

void print_mat(double* p_mat, ///<matrix to be printed, column-major
					  int n_row, 
					  int n_col, 
					  char* name) {
	if (p_mat == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s matrix (double-precision floats): (%u x %u)\n", 
		name, n_row, n_col);
	for (int ir = 0; ir < n_row; ++ir) {
		for (int ic = 0; ic < n_col; ++ic) {
			printf("%4.2f ", p_mat[ir + n_row * ic]);
		}
		printf("\n");
	}
}

static void print_lmat(double* p_mat, ///<matrix to be printed, column-major
					  int n_row, 
					  int n_col, 
					  char* name) {
	// print lower triangular matrix
	if (p_mat == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s matrix (double-precision floats): (%u x %u)\n", 
		name, n_row, n_col);
	for (int ir = 0; ir < n_row; ++ir) {
		for (int ic = 0; ic <= ir; ++ic) {
			printf("%4.2f ", p_mat[ir + n_row * ic]);
		}
		printf("\n");
	}
}

static void print_vec(int* p_vec, ///<vector to be printed
					  int n, ///< number of elements 
					  char* name) {
	if (p_vec == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s vector (positive integers): (n = %u)\n", 
		name, n);
	for (int i = 0; i < n; ++i) {
		printf("%4u ", p_vec[i]);
	}
	printf("\n");
}
static void print_vec(double* p_vec, ///<vector to be printed
					  int n, ///< number of elements 
					  char* name) {
	if (p_vec == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s vector (double-precision floats): (n = %u)\n", 
		name, n);
	for (int i = 0; i < n; ++i) {
		printf("%4.2f ", p_vec[i]);
	}
	printf("\n");
}



void VarModel_dump(VarModel* varm) {
	if (varm == NULL) {
		printf("VarModel object not created.\n");
		return;
	}
	printf("VarModel object %p, %u dimensions, "
		"%u seasonal periods, lookback = %d.\n",
		varm, varm->n, varm->ns, varm->lb);

	print_vec(varm->mu, varm->n, "mu");
	//print_mat(varm->phi, varm->n, varm->n*varm->phi_cnt, "phi");
	int n = varm->n;
	for (int i = 0; i<varm->npos; ++i) {
		print_mat(&varm->poly[i*n*n], n, n, "poly");
	}
	print_vec(varm->pos, varm->npos, "pos");

	print_mat(varm->cov, varm->n, varm->n, "cov");
	print_lmat(varm->lcov, varm->n, varm->n, "lcov");

	printf("Determinant of the covariance matrix: %g\n", 
		SQR(varm->sqrt_d));

}

VarModel_Err VarModel_setCov(VarModel* varm, 
							 double* cov_in, int cov_len) {
	int dim = varm->n;
	if (cov_len != dim * dim) {
		return VARMODEL_COV_MATRIX_ERROR;
	}

	//printf("Setting cov...\n");
	memcpy(varm->cov, cov_in, cov_len*sizeof(double));
	memcpy(varm->lcov, cov_in, cov_len*sizeof(double));
	//print_mat(varm->cov, dim, dim, "cov");


	integer n = dim;
	integer info;
	//printf("Setting lcov...\n");
	dpotrf_("L", &n, varm->lcov, &n, &info);

	if (info == 0) {
		//print_lmat(varm->lcov, dim, dim, "lcov");
		varm->sqrt_d = 1.0;
		for (int i=0; i<dim; ++i) {
			varm->sqrt_d *= varm->lcov[i + dim*i];
		}
		//printf("Sqrt(determinant) = %8.4f\n", varm->sqrt_d);
	} else {
		if (info < 0) {
			//printf("%d-th parameter of Cholesky Fact. is illegal\n", -info);
		} else  {
			//printf("The cov matrix is not SPD\n");
		}
		return VARMODEL_COV_MATRIX_ERROR;
	}

	return VARMODEL_OK;
}

VarModel_Err VarModel_estMu(VarModel* varm, double* panel, int panel_dim,
							int panel_len) {
	if (panel_dim != varm->n) {
		printf("Error dimensions of the panel data.\n");
		return VARMODEL_DIM_ERR;
	}

	//compute the average
	integer n = panel_dim;
	double a = 1;
	integer inc = 1;
	for (int i=0; i<panel_len; ++i) { //sum
		daxpy_(&n, &a, &panel[i*panel_dim], &inc, varm->mu,  &inc);
	}

	for (int i=0; i<panel_dim; ++i) {
		varm->mu[i] /= panel_len;
	}

	return VARMODEL_OK;
}



VarModel_Err VarModel_logLikelihood(VarModel* varm, double* panel,
			int panel_dim, int panel_len, double* logl_out, 
			double* a_out, int n_a) {

	int nea = panel_len - varm->lb;
	if (nea <= 0) {
		//printf("Not enough samples in the panel data.\n");
		return VARMODEL_NOT_ENOUGH_SAMPLE;
	}

	if (a_out != NULL && n_a != nea) {
		return VARMODEL_ARRAY_ERROR;
	}

	double tl; // temp variable for log-likelihood
	tl = nea * ( -0.5*panel_dim*log(2*_PI) - log(varm->sqrt_d));

	// two temporary arrays
	double* w = (double*)calloc(panel_dim, sizeof(double));
	double* a = (double*)calloc(panel_dim, sizeof(double));

	integer n = panel_dim;
	integer inc = 1;
	double alpha, beta;
	int ia_out = 0;

	for (int t=varm->lb; t<panel_len; ++t) {
		//log-likelihood for each realization of white noise
		memset(a, 0, panel_dim*sizeof(double));

		for (int ip=0; ip<varm->npos; ++ip) {
			// w_t = y_t - mu
			memcpy(w, &panel[panel_dim * (t - varm->pos[ip])],
				   panel_dim * sizeof(double));
			alpha = -1;
			daxpy_(&n, &alpha, varm->mu, &inc, w, &inc);

			//compute with char. polyn. 
			alpha = 1; beta = 1;
			//dgemv_("n", &n, &n, &alpha, &varm->poly[n*n*varm->pos[ip]], &n,
			dgemv_("n", &n, &n, &alpha, &varm->poly[n*n*ip], &n,
				   w, &inc, &beta, a, &inc);
		}
		// output a
		if (a_out != NULL) {
			memcpy(&a_out[(ia_out++)*panel_dim], a, sizeof(double)*panel_dim);
		}
		
		// from pdf of multi-normal dist.
		// compute (L^T)^-1 * a
		integer nrhs = 1;
		integer info;
		dtrtrs_("L", "N", "N", &n, &nrhs, varm->lcov, &n, a, &n, &info);
		double nm;
		nm = dnrm2_(&n, a, &inc);
		tl += -0.5*nm*nm;
	}

	(*logl_out) = tl;
	free(w);
	free(a);

	return VARMODEL_OK;
}


VarModel_Err VarModel_estimate(VarModel* vm, double* panel, int dim, int len, double* diff, double* mdll) {

	//double *tp = (double*)malloc(200*sizeof(double));
	*diff = 0;

	if (dim != vm->n) return VARMODEL_DIM_ERR;
	if (len <= vm->lb) return VARMODEL_NOT_ENOUGH_SAMPLE;
	if (vm->npos <= 0) return VARMODEL_NO_AR_PARA;

	// update mean
	//memset(vm->mu, 0, sizeof(double)*dim);
	double mu0;
	for (int i=0; i<dim; ++i) {
		mu0 = vm->mu[i];
		vm->mu[i] = 0;
		for (int t=0; t<len; ++t) {
			vm->mu[i] += panel[i + t*dim];
		}
		vm->mu[i] /= len;
		*diff += SQR(mu0-vm->mu[i]);
	}

	// construct matrix A in Ax=b regression
	int T = len - vm->lb;
	int p = vm->npos - 1;
	int m = T*dim;  // no. of rows;
	int n = dim*dim*p;  // no. of columns
	if (m<n) return VARMODEL_NOT_ENOUGH_SAMPLE;

	double *A = (double*)calloc(m*n, sizeof(double));
	double* b = (double*)calloc(dim*T, sizeof(double));

	for (int j = 0; j < T; ++j) {
		for (int ipos = 1; ipos < vm->npos; ++ipos) {
			for (int idim = 0; idim < dim; ++idim) {
				//element of A
				int i = (ipos-1)*dim + idim;
				double a = panel[(vm->lb - vm->pos[ipos] + j)*dim + idim] - vm->mu[idim]; 
				for (int k=0; k<dim; ++k) {
					A[j*dim+k + T*dim * (i*dim + k)] = a; //A' (*) I
				}
			}
		}
	}

	// adjust panel by mean values
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < T; ++j) {
			b[i + j*dim] = panel[i + (vm->lb+j)*dim] - vm->mu[i];
		}
	}


	// LS solver: lapack dgelss()
	integer m1 = m;
	integer n1 = n;
	integer nrhs = 1;
	double rcond = 1e5;  //
	integer rank = 0;
	double* s = (double*)calloc(n, sizeof(double)); //singular values
	integer info;

	// workspace query first, check mkl. man.
	integer lwork = -1; 
	double* work = (double*)calloc(1, sizeof(double));;  //temp array

	dgels_("N", &m1, &n1, &nrhs, A, &m1, b, &m1, work, &lwork, &info);
	//dgelss_(&m1, &n1, &nrhs, A, &m1, b, &m1, s, &rcond, &rank, work, &lwork, &info);  

	// LS regression
	lwork = integer(work[0]);
	work = (double*)realloc(work, sizeof(double)*lwork);
	dgels_("N", &m1, &n1, &nrhs, A, &m1, b, &m1, work, &lwork, &info);
	//dgelss_(&m1, &n1, &nrhs, A, &m1, b, &m1, s, &rcond, &rank, work, &lwork, &info);  

	// update vm->poly
	for (int i = 0; i< n; ++i ) {
		*diff += SQR(b[i]-vm->poly[dim*dim+i]);
	}
	memcpy(&vm->poly[dim*dim], b, sizeof(double)*n);
	
	// update cov and lcov
	double dy_ll;
	double* rsd = (double*)calloc(T*dim, sizeof(double));
	VarModel_logLikelihood(vm, panel, dim, len, &dy_ll, rsd, T);
	double* cov = (double*)calloc(dim*dim, sizeof(double));
	for (int i = 0; i<dim; ++i) {
		for (int j=0; j<dim; ++j) {
			for (int ir = 0; ir < T; ++ir) {
				cov[i + j*dim] += rsd[i + dim*ir] * rsd[j + dim*ir];
			}
			cov[i + j*dim] /= T;
		}
	}

	for (int i=0; i<dim*dim; ++i) {
		*diff += SQR(vm->cov[i]-cov[i]);
	}
	VarModel_setCov(vm, cov, dim*dim);
	// compute the final max dmll
	VarModel_logLikelihood(vm, panel, dim, len, mdll);

	free(work);
	free(b);
	free(A);
	return VARMODEL_OK;
}