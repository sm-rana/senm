#include "limits.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "SenMCoreIncs.h"
#include "clapack_3.2.1\clapack.h"
#include "VarModel.h"

#include "levmar.h" // Levenberg-Marquardt algorithm ( masud-- 6/23/2016)

//VarModel_new(dim=7, ns=1, s=1, p=[2,1], &vm);
VarModel_Err VarModel_new(int n_in, int ns_in, int* s_in, int* p_in, VarModel** varm_out) 
{
	VarModel* p = (VarModel*)calloc(1, sizeof(VarModel)); // p = varm
	if (p == NULL) return VARMODEL_ALLOC_ERROR;

	p->n = n_in; //=7
	p->ns = ns_in; //=1
	p->s = (int*)calloc(ns_in + 1, sizeof(int)); //=2
	p->p = (int*)calloc(ns_in + 1, sizeof(int)); //=2
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

//VarModel_setPhi(vm, phi, dim*dim*np);
VarModel_Err VarModel_setPhi(VarModel* varm, 
							 double* phi_in, 
							 int phi_size) {
	long n = (long)varm->n; //=7
	int phi_cnt = varm->phi_cnt;
	if (phi_size != n * n * (phi_cnt)) 
		return VARMODEL_PHI_ARRAY_ERROR;
	memcpy(varm->phi, phi_in, sizeof(double)*phi_size);
	
	double* tp_poly; // tp_poly[245] temporary array of matrices
	int* tp_pos; // masud tp_pos[6]// temporary array of non-zero para positions
	int n_items = 1; //= 6;
	for (int i_ns = 0; i_ns <= varm->ns; ++ i_ns) {
		n_items *= (varm->p[i_ns] + 1);
	}

	tp_poly = (double*)calloc(n*n*n_items, sizeof(double)); 
	tp_pos = (int*)calloc(n_items, sizeof(int)); //tp_pos[6]
	varm->poly = (double*)calloc(n*n*n_items, sizeof(double)); 
	varm->pos = (int*)calloc(n_items, sizeof(int)); //
	for (int i = 0; i < n; ++i) {
		tp_poly[i + n*i] = 1; // tp_poly[0,8,16,... 49] // diagonals are made 1
	} // iden. mat 
	int len = 1; // =1; length of the tp_pos (also no. matrices in tp_poly)

	double alpha = -1; //= -1
	double beta = 0; // =0
	int ip0 = 0; //=0
	for (int i_ns = 0; i_ns <= varm->ns; ++i_ns) {
		int p = varm->p[i_ns];//=1,2
		int s = varm->s[i_ns];//=24,25?

		for (int ip =0; ip < p; ++ip) { // new polynomial
			for (int itp = 0; itp < len; ++itp) { 
				// matrix multiplication
				//dgemm(): tp_poly = -1*tp_poly*phi
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
							&varm->poly[p_poly*n2], &inc); //poly[]=poly[] + tp_poly[]
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

//
//masud: print_mat(buf, dim, buf_size, "Buffer");
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
	// n_row = 7
	// n_col = 194
	// p_mat = buf

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

void VarModel_dump2(VarModel* varm, int iteration, FILE ** filePointer) {
	if (varm == NULL) {
		printf("VarModel object not created.\n");
		return;
	}

	FILE *fp; //masud: debug_varm.txt
	fp = *filePointer;

	int n = varm->n; //=7

	fprintf(fp,"\n--------------------------\nIteration = %d\n------------\n", iteration);
	fprintf(fp,"VarModel object %p\n %u dimensions, "
		"%u seasonal periods, lookback = %d.\n",
		varm, varm->n, varm->ns, varm->lb);

	//print_vec(varm->mu, varm->n, "mu");
	//print_mat(varm->phi, varm->n, varm->n*varm->phi_cnt, "phi");

	fprintf(fp,"mu:\t");
	for (int i = 0; i < n; i++)
		fprintf(fp, "%.3f\t", varm->mu[i]);
	fprintf(fp, "\n\n");
		
	for (int i = 0; i<varm->npos; ++i) 
	{
		//print_mat(&varm->poly[i*n*n], n, n, "poly");
		fprintf(fp, "\nPhi Matrix %d:\n", i);
		for (int ir = 0; ir <n; ++ir)
		{
			for (int ic = 0; ic < n; ++ic) {
				fprintf(fp, "%.4f  ", varm->poly[i*n*n + (ir + n * ic)]);
			}
			fprintf(fp,"\n");
		}
	}

	fprintf(fp, "\nCovariance Matrix:\n");
	print_mat(varm->cov, varm->n, varm->n, "cov");
	for (int ir = 0; ir <n; ++ir)
	{
		for (int ic = 0; ic < n; ++ic) {
			fprintf(fp, "%.4f  ", varm->cov[ir + n * ic]);
		}
		fprintf(fp, "\n");
	}

	//print_lmat(varm->lcov, varm->n, varm->n, "lcov");

	//printf("Determinant of the covariance matrix: %g\n",
		//SQR(varm->sqrt_d));

}

VarModel_Err VarModel_setCov(VarModel* varm, 
							 double* cov_in, int cov_len) {
	int dim = varm->n;
	if (cov_len != dim * dim) {
		return VARMODEL_COV_MATRIX_ERROR;
	}

	//printf("Setting cov...\n");
	memcpy(varm->cov, cov_in, cov_len*sizeof(double)); //lcov[7x7]
	memcpy(varm->lcov, cov_in, cov_len*sizeof(double)); //lcov[7x7]
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
	
	//masud: panel_len is 1 in senm.cpp
	//masud: VarModel_estMu(vm, mu, dim, 1);
	if (panel_dim != varm->n) {
		printf("Error dimensions of the panel data.\n");
		return VARMODEL_DIM_ERR;
	}

	//compute the average
	integer n = panel_dim; //=7
	double a = 1;
	integer inc = 1;
	for (int i=0; i<panel_len; ++i) { //sum
		daxpy_(&n, &a, &panel[i*panel_dim], &inc, varm->mu,  &inc); //mu = mu + panel[] // mu = varm->mu
	}
	

	for (int i=0; i<panel_dim; ++i) {
		varm->mu[i] /= panel_len;
	}

	return VARMODEL_OK;
}


//VarModel_logLikelihood(1 vm, 2 &buf[t*dim], 3 dim, 4 vm->lb+1, 5 &dm_ll);		//from main()
//VarModel_logLikelihood(vm, panel, dim, len, &dy_ll, rsd, T_window);			//from VarModel_estimate()
VarModel_Err VarModel_logLikelihood(VarModel* varm, double* panel,
			int panel_dim, int panel_len, double* logl_out, 
			double* a_out, int n_a) {
	// masud:
	//1 varm: object of class VarModel
	//2 panel -- pointer to buf -- current demand 
	//3 panel_dim = dim = 7 
	//4 panel_len = lookback +1 = 27
	//5 dm_ll = model likelyhood
	// n_a = 0; set in the header file and has the description: /*effective no. w. noise*/

	int nea = panel_len - varm->lb; //masud: nea = 1
	if (nea <= 0) {
		//printf("Not enough samples in the panel data.\n");
		return VARMODEL_NOT_ENOUGH_SAMPLE;
	}

	if (a_out != NULL && n_a != nea) {
		return VARMODEL_ARRAY_ERROR;
	}
	
	//masud: panel_dim = 7, nea = 1
	///sqrt_d < sqrt of cov determinant
	double tl; // temp variable for log-likelihood
	tl = nea * (-0.5*panel_dim*log(2*_PI) - log(varm->sqrt_d));
	

	// two temporary arrays
	double* w = (double*)calloc(panel_dim, sizeof(double)); //masud: w[7] 
	double* a = (double*)calloc(panel_dim, sizeof(double)); // a[7]

	integer n = panel_dim; // = 7
	integer inc = 1;	// = 1
	double alpha, beta;
	int ia_out = 0;	

	//for(t = 26; t<27; t++}
	for (int t=varm->lb; t<panel_len; ++t) 
	{
		//log-likelihood for each realization of white noise
		memset(a, 0, panel_dim * sizeof(double)); //a[0:6] = 0;

		for (int ip=0; ip<varm->npos; ++ip) 
		{ // npos = 6; pos[] = 0,1,2,24,25,26
			
			// masud: w contains demand at a certain lag
			memcpy(w, &panel[panel_dim * (t - varm->pos[ip])],
				   panel_dim * sizeof(double)); //masud: panel=buf[7*(26-pos[0:npos])] 
			
			alpha = -1;
			daxpy_(&n, &alpha, varm->mu, &inc, w, &inc); //paulo: w = w - mu; (w is nothing but demand estimates minus the mean)
			//daxpy(N,DA,DX,INCX,DY,INCY) == dy(iy) = dy(iy) + da*dx(ix)

			//compute with char. polyn. 
			alpha = 1; beta = 1;
			//dgemv_("n", &n, &n, &alpha, &varm->poly[n*n*varm->pos[ip]], &n,

			//paulo: a = poly[n*n*ip]*w + a // x_t = phi*x_(t-1) + a
			dgemv_("n", &n, &n, &alpha, &varm->poly[n*n*ip], &n,
				   w, &inc, &beta, a, &inc);
									
		}
		// output a
		
		//masud: the following doesn't run
		if (a_out != NULL) {
			memcpy(&a_out[(ia_out++)*panel_dim], a, sizeof(double)*panel_dim);
		}
		

		// from pdf of multi-normal dist.
		// compute (L^T)^-1 * a
		integer nrhs = 1; //=1
		integer info;
		dtrtrs_("L", "N", "N", &n, &nrhs, varm->lcov, &n, a, &n, &info);
		double nm;
		nm = dnrm2_(&n, a, &inc);
		tl += -0.5*nm*nm;	//masud: time series likelihood
	}

	(*logl_out) = tl;
	free(w);
	free(a);

	return VARMODEL_OK;
}

//masud: call from senm -- VarModel_estimate(vm, buf, dim, buf_size, &diff, &mdll);
// buf_size = 26 + 168 = 194 == len
VarModel_Err VarModel_estimate(VarModel* vm, double* panel, int dim, int len, double* diff, double* mdll) {

	//double *tp = (double*)malloc(200*sizeof(double));
	*diff = 0;

	if (dim != vm->n) return VARMODEL_DIM_ERR;
	if (len <= vm->lb) return VARMODEL_NOT_ENOUGH_SAMPLE;	// = 194
	if (vm->npos <= 0) return VARMODEL_NO_AR_PARA;

	// update mean
	//memset(vm->mu, 0, sizeof(double)*dim);
	//masud: panel == buf
	double mu0;
	for (int i=0; i<dim; ++i) {
		mu0 = vm->mu[i];
		vm->mu[i] = 0;
		for (int t=0; t<len; ++t) { //len = 194
			vm->mu[i] += panel[i + t*dim];
		}
		vm->mu[i] /= len;
		*diff += SQR(mu0-vm->mu[i]);
	}

	// construct matrix A in Ax=b regression
	int T = len - vm->lb; //masud: T = 194-26 = 168
	int p = vm->npos - 1; // p = 5 == number of phi matrices
	int m = T*dim;  //masud: no. of rows; m[168][7]
	int n = dim*dim*p;  //masud: = 49*5 total number of elements in the phi matrices... 
	if (m<n) return VARMODEL_NOT_ENOUGH_SAMPLE;

	double *A = (double*)calloc(m*n, sizeof(double));	//A[168][5] = lb, n = p -- p = # of AR params
	double* b = (double*)calloc(dim*T, sizeof(double));	// b[168]  
		
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

	//debug: masud -- 
	for (int jj = 0; jj < 3*49; jj++)
	{
		double temp = -1 * b[jj];
		b[jj] = temp;
	}

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


//masud: call from senm -- VarModel_estimate(vm, buf, dim, buf_size, &diff, &mdll);

void Var_Function(double *phi, double *y, int m1, int n, void *data)
{	/*
	WARNING: BE CAREFUL ABOUT USING THIS FUNCTION
	//masud: 6/23/2016
	// This is the model that the Levenberg-Marquardt function uses to fit the parameters. 
	
	####### WARNING: This function needs to be changed if the Time-Series Model structure is changed ##########
	// phi	-->	vector containing initial guess of the parameters
	// y	-->	model estimate
	// m1	-->	number of parameters to be estimated !! NOTE: this number cannot be used for matrix operation because we have additional product terms
	// n	-->	number of data points
	// data -->	this is any data that is needed by the function. I am passing the  past observation of the time series and total phis which includes product terms (=5 for now). 
	//			The	type of data has to be "void" and hence need to be typcasted to the appropriate type before being used. 
	**** NOTE ****
	// This function is to be passed to dlevmar_dif() or dlevmar_der() function found in levmar.lib 
	// For more details on the Levenberg-Marquardt algorithm visit:  http://users.ics.forth.gr/~lourakis/levmar/ 
	*/
	double ** data_mat = (double **)data;		// past time-series observation
	double m = data_mat[0][0]; // total number of phi + product terms = 20
	double* A_matrix = data_mat[1];	// n row x m column matrix		//e.g., n = 142*dim, m = 5*dim

	int dim = data_mat[0][2];
	int np = data_mat[0][3]; //np = 3

	integer nrow = (integer)n;
	integer mcol = (integer)m;
	integer incx = 1; // increment of X and Y vector 
	
	//todo 5: the product phi terms need to be automatically determined using p and s values;
	double* phi1 = (double*)calloc(m, sizeof(double)); // phi1[20] 

	for (int ii = 0; ii < dim*np; ii++)
		phi1[ii] = phi[ii];

	for (int ii = 0; ii < dim; ii++)
	{
		phi1[np*dim + ii] = - phi[dim *0 + ii] * phi[(np - 1)*dim + ii];
		phi1[(np+1)*dim + ii] = - phi[dim * 1 + ii] * phi[(np - 1)*dim + ii];
	}

	/*FILE* fp_var_func;
	fopen_s(&fp_var_func, "debug_Var_Func.txt", "w");

	fprintf(fp_var_func, "Phi Var_Function dgemv_():\n");
	for (int ii = 0; ii < 5; ii++)
	{
		for (int jj = 0; jj < 4; jj++)
			fprintf(fp_var_func, "%f ", phi1[ii * 4 + jj]);
		fprintf(fp_var_func, "\n");
	}*/

	double alpha = 1;
	double beta = 0;
	dgemv_("N", &nrow, &mcol, &alpha, A_matrix, &nrow, phi1, &incx, &beta, y, &incx);
	
	/*fprintf(fp_var_func, "\n\nY-estimates:\n");
	for (int ii = 0; ii < 46; ii++)
	{
		for (int jj = 0; jj < 4; jj++)
			fprintf(fp_var_func, "%f ", y[ii*4 + jj]);
		fprintf(fp_var_func, "\n");
	}
	fprintf(fp_var_func, "\n\n");
	fclose(fp_var_func);*/
	free(phi1);
}


//panel = buf
//dim = 7
//len = 194
VarModel_Err VarModel_estimate2(VarModel* vm, double* panel, int dim, int len, double* diff, double* mdll) 
{
	*diff = 0;
	if (dim != vm->n) return VARMODEL_DIM_ERR;
	if (len <= vm->lb) return VARMODEL_NOT_ENOUGH_SAMPLE;	// = 194
	if (vm->npos <= 0) return VARMODEL_NO_AR_PARA;

	// update mean
	//memset(vm->mu, 0, sizeof(double)*dim);
	//masud: panel == buf
	double mu0;
	for (int i = 0; i<dim; ++i) {
		mu0 = vm->mu[i];
		vm->mu[i] = 0;
		//masud: 6-16-2016: changed the averaging from t = 0 onward. Jinduan used to do 
		//     it from t = -26 onward where -26th data is just last 26 data copied to the back
		for (int t = vm->lb; t<len; ++t) { //len = 194
			vm->mu[i] += panel[i + t*dim];
		}
		vm->mu[i] /= (len-vm->lb);
		*diff += SQR(mu0 - vm->mu[i]);
	}
	/*
	 * masud:
	 * todo 4: clean up VarModel_estimate2() --VarModel.cpp. Too many vectors declared. Not memory efficient. 
	 * VarModel estimate being the least expensive part -- who cares of inefficiency!!
	 * construct matrix A in Ax=b regression
	 * dim == the unknown demand variables (i.e., number of clusters for grouped nodes)
	 */
	

	int T_window = len - vm->lb;				//masud: T = 194-26 = 168
	int total_phi = vm->npos - 1;				// p = 5 == number of phi matrices
	int rowA = T_window*dim;					//m = [168*7]; masud: no. of rows; m[168][7]
	int rowA_uni = (T_window - vm->lb)*dim;		//m = [142*7]; masud: no. of rows; m[142][7]
	//masud: 6-12-2016. Why did I exclude first 26 hrs? -- (6-16-2016)Because we can't estimate parameter using those values... you idiot.
	int n = dim*dim*total_phi;					//n = [5*7*7] masud: = total number of elements in the phi matrices... 
	int n_univariate = total_phi * dim * 1;		// = 5 * 7 * 1
	if (rowA<n) return VARMODEL_NOT_ENOUGH_SAMPLE;

	double *A = (double*)calloc(rowA*n, sizeof(double));	//A[168*7][5*7*7] = lb, n = p -- p = # of AR params
	double *A_univariate = (double*)calloc(rowA_uni*n_univariate, sizeof(double)); //masud: for univariate case... ignoring the off-diagonals
	double *A_uni_lm = (double*)calloc(rowA_uni*n_univariate, sizeof(double)); //masud: for univariate case... ignoring the off-diagonals

	double* b_multi = (double*)calloc(dim*T_window, sizeof(double));	// b[168*7]  
	double* b_uni = (double*)calloc(rowA_uni, sizeof(double));	// b_uni[142*7]  
	double* b_uni_lm = (double*)calloc(rowA_uni, sizeof(double));	// b_uni[142*7] -- lm = Levenberg - Marquardt algorithm

	double* b = (double*)calloc(dim*T_window, sizeof(double));	// b[168*7]  

	FILE *fp; //"debug_EstimatePhi.txt"
	fopen_s(&fp, "Debug_EstimatePhi.txt", "w");
	
	for (int j = 0; j < T_window; ++j) 
	{
		for (int ipos = 1; ipos < vm->npos; ++ipos) 
		{
			for (int idim = 0; idim < dim; ++idim) 
			{
				//element of A
				int i = (ipos - 1)*dim + idim;
				double a = panel[(vm->lb - vm->pos[ipos] + j)*dim + idim] - vm->mu[idim];
				for (int k = 0; k<dim; ++k) 
				{
					A[j*dim + k + T_window*dim * (i*dim + k)] = a; //A' (*) I
				}
			}
		}
	}
	
	int skip = vm->lb * dim; //skip = 182; //cluster = 78
	int countt = 0;
	for (int nphi = 0; nphi < total_phi; nphi++)
	{
		for (int column = 0; column < dim; column++)
		{
			int hopscotch = (nphi * dim*dim * rowA) + (rowA * (dim + 1) * column);
			for (int i = skip; i < rowA; i++)
			{
				A_uni_lm[countt] = A[i + hopscotch];
				A_univariate[countt] = A[i + hopscotch];
				countt++;
			}
		}
	}

	//fprintf(fp, "\n\n\nb vector:\n");
	// adjust panel by mean values
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < T_window; ++j) {
			b_multi[i + j*dim] = panel[i + (vm->lb + j)*dim] - vm->mu[i];
		}
	} 

	countt = 0;
	for (int i = skip; i < rowA; i++)
	{
		b_uni_lm[countt] = b_multi[i];
		b_uni[countt] = b_multi[i];
		countt++;
	}

	fprintf(fp, "A Matrix: \n");
	for (int i = 0; i < rowA_uni; i++)
	{
		for (int j=0; j<dim * total_phi; j++)
		{
			fprintf(fp, "%f ", A_univariate[j * rowA_uni + i]);
		}
		fprintf(fp, "\n");
	}

	countt = 0;
	fprintf(fp, "\n\nb Matrix:\n");
	for (int i = 0; i < rowA_uni; i++)
		fprintf(fp, "%f\n", b_uni[countt++]);

	// LS solver: lapack dgelss()
	integer m1 = rowA; //m1[168*7] -- rows of A
	integer m1_uni = rowA_uni; //m1[168*7] -- rows of A
	integer n1 = n; // n = [5*7*7] -- columns of A
	integer n1_uni = n_univariate; // n = [5*7*1] -- columns of A
	integer nrhs = 1; // number of right hand side
	double rcond = 1e5;  //
	integer rank = 0;
	double* s = (double*)calloc(n, sizeof(double)); //singular values
	integer info;

	/*
		The following is a linear approaximation of time sries model parameters. Which means that the product terms are considered as different parameters
		instead of actual combination of parameters. For example, phi1*phi2 is considered a unique parameter. 
	*/

	integer lwork = -1;
	double* work_uni = (double*)calloc(1, sizeof(double));;  //temp array
	dgels_("N", &m1_uni, &n1_uni, &nrhs, A_univariate, &m1_uni, b_uni, &m1_uni, work_uni, &lwork, &info);
	if (info != 0)
		printf("\nProblem in linear Estimation of Phi\n");
	// LS regression
	lwork = integer(work_uni[0]);
	work_uni = (double*)realloc(work_uni, sizeof(double)*lwork);
	dgels_("N", &m1_uni, &n1_uni, &nrhs, A_univariate, &m1_uni, b_uni, &m1_uni, work_uni, &lwork, &info);
	if (info != 0)
		printf("\nProblem in linear Estimation of Phi\n");
	/* 
	////%%%%%%%%%%%%%%%%%%%%%% NON-LINEAR ESTIMATION OF PHI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Masud: Non-linear estimation of phi using linear estimation as an initial guess:
		- I need to send A_univariate and b_uni to levmar algorithm:
		- b_uni_lm == B vector
		- phi_lm == phi vector (dgels_() -- inserted linear estimate of phi in b_uni)
		- A_uni_lm is is the A matrix
	*/

	double *data[2];
	const int temp = 4;
	double data1[temp];
	data1[0] = (double)n1_uni;	// columns
	data1[1] = (double)m1_uni;	// rows
	data1[2] = dim;	// dim = 4 (unknown demand to be estimated)
	data1[3] = (double) vm->phi_cnt; //

	data[0] = data1;
	data[1] = A_uni_lm;
	int np = vm->phi_cnt; //np = 3

	void(*pVar_func)(double*, double*, int, int, void*); // pointer to Var_func
	pVar_func = &Var_Function;

	// masud: column and row symbols have changed here -- don't bang your head.
	// before n = col, m = row, for levmar n = row, m = col
	int m_params;
	m_params = np * dim; //=12
	
	int n_data = rowA_uni;
	int max_trials = 1000;

	double opts[LM_OPTS_SZ], infor[LM_INFO_SZ]; //opts[5]; info[10];
	opts[0] = LM_INIT_MU; 
	opts[1] = 1E-15; 	opts[2] = 1E-15; 	opts[3] = 1E-20;
	opts[4] = -1; // masud: Jacobial will be calculated using central difference
	
	fprintf(fp, "\n\nb_uni -- linear estimate:\n");
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			fprintf(fp, "%f ", b_uni[i*dim + j]);
			//b_uni[i*dim + j] = b_uni[i*dim + j] - 0.5;
		}
		fprintf(fp, "\n");
	}

	double* phi_lm = (double*)calloc(np*dim, sizeof(double));

	fprintf(fp, "\n\nphi_lm Initial Guess:\n");
	for (int i = 0; i < np; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			phi_lm[i*dim + j] = b_uni[i*dim + j];
			fprintf(fp, "%f ", phi_lm[i*dim + j]);
		}
		fprintf(fp, "\n");
	}
				  
	// b_uni == phi vector 
	dlevmar_dif(pVar_func, phi_lm, b_uni_lm, m_params, n_data, max_trials, opts, infor, NULL, NULL, data); // without Jacobian
	printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", infor[5], infor[6], infor[1], infor[0]);

	printf("Info Vector:\n");
	for (int i = 0; i < 10; i++)
	{
		printf("[%d] = %f\n",i, infor[i]);
	}

	fprintf(fp, "\n\nphi_lm AFTER NON-LIN ESTIMATION:\n");
	for (int i = 0; i < np; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			fprintf(fp, "%f ", phi_lm[i*dim + j]);
			b_uni[i*dim + j] = phi_lm[i*dim + j];
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	/*countt = 0;
	fprintf(fp, "\n\nb_uni_lm Matrix:\n");
	for (int i = 0; i < rowA_uni; i++)
		fprintf(fp, "%f %f\n", b_uni_lm[countt],b_uni[countt++]);
	fclose(fp);*/

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// calculating the product terms of phi
	for (int ii = 0; ii < dim; ii++)
	{
		b_uni[np*dim + ii] = -phi_lm[dim * 0 + ii] * phi_lm[(np - 1)*dim + ii];
		b_uni[(np + 1)*dim + ii] = -phi_lm[dim * 1 + ii] * phi_lm[(np - 1)*dim + ii];
	}
	free(phi_lm);

	for (int i = 0; i < dim*T_window; i++)
		b[i] = 0.0;

	countt = 0;
	for (int nphi = 0; nphi < total_phi; nphi++)
	{
		for (int column = 0; column < dim; column++)
		{
			int hopscotch = (nphi * dim*dim) + ((dim+1) * column);
			b[hopscotch] = -1 * b_uni[countt++];
		}
	}
	
	for (int i = 0; i< n; ++i) {
		*diff += SQR(b[i] - vm->poly[dim*dim + i]); 
	}
	memcpy(&vm->poly[dim*dim], b, sizeof(double)*n);

	double dy_ll;
	double* rsd = (double*)calloc(T_window*dim, sizeof(double)); // rsd [168*7]
	VarModel_logLikelihood(vm, panel, dim, len, &dy_ll, rsd, T_window);
	double* cov = (double*)calloc(dim*dim, sizeof(double));
	
	for (int i = 0; i<dim; ++i) 
	{
		for (int j = 0; j<dim; ++j) 
		{
			if (i != j)	//masud: no need to calculate covariances but only variance
				continue;	//
			for (int ir = 0; ir < T_window; ++ir) 
			{
				cov[i + j*dim] += rsd[i + dim*ir] * rsd[j + dim*ir];
			}
			cov[i + j*dim] /= T_window;
		}
	}

	for (int i = 0; i<dim*dim; ++i) {
		*diff += SQR(vm->cov[i] - cov[i]);
	}
	VarModel_setCov(vm, cov, dim*dim);
	// compute the final max dmll
	VarModel_logLikelihood(vm, panel, dim, len, mdll);

	free(work_uni);
	//free(work);
	free(b);
	free(b_uni);
	free(b_multi);

	free(A);
	free(A_univariate);
	free(A_uni_lm);
	return VARMODEL_OK;
}



int VarModel_estimate_levmar(void(*Var_func)(double *p, double *hx, int m, int n, void *adata),
	double*xt, double*xt_1, double* params0, int m_params, int n_data, int max_trials, double* opts, double* info)
{
	//Masud : is this a garbage function? Do we need this? (11/04/2016)
	int ret;
	// optimization control parameters; passing to levmar NULL instead of opts reverts to defaults 

	//ret = dlevmar_der(expfunc, jacexpfunc, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // with analytic Jacobian
	ret = dlevmar_dif(Var_func, params0, xt, m_params, n_data, max_trials, opts, info, NULL, NULL, xt_1); // without Jacobian

	printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
	printf("Best fit parameters: %.7g %.7g %.7g\n", params0[0], params0[1], params0[2]);

	return ret;

	/*
	* Similar to dlevmar_der() except that the Jacobian is approximated internally with the aid of finite differences.
	* Broyden's rank one updates are used to compute secant approximations to the Jacobian, effectively avoiding to call
	* func several times for computing the finite difference approximations.
	* If the analytic Jacobian is available, use dlevmar_der() above.
	*
	* Returns the number of iterations (>=0) if successful, -1 if failed
	*
	*/

	//int dlevmar_dif(
	//	void(*func)(double *p, double *hx, int m, int n, void *adata), /* functional relation describing measurements.
	//																   * A p \in R^m yields a \hat{x} \in  R^n
	//																   */
	//	double *p,         /* I/O: initial parameter estimates. On output contains the estimated solution */
	//	double *x,         /* I: measurement vector. NULL implies a zero vector */
	//	int m,             /* I: parameter vector dimension (i.e. #unknowns) */
	//	int n,             /* I: measurement vector dimension */
	//	int itmax,         /* I: maximum number of iterations */
	//	double opts[5],    /* I: opts[0-4] = minim. options [\tau, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
	//					   * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and the
	//					   * step used in difference approximation to the Jacobian. If \delta<0, the Jacobian is approximated
	//					   * with central differences which are more accurate (but slower!) compared to the forward differences
	//					   * employed by default. Set to NULL for defaults to be used.
	//					   */
	//	double info[LM_INFO_SZ],
	//	/* O: information regarding the minimization. Set to NULL if don't care
	//	* info[0]= ||e||_2 at initial p.
	//	* info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, \mu/max[J^T J]_ii ], all computed at estimated p.
	//	* info[5]= # iterations,
	//	* info[6]=reason for terminating: 1 - stopped by small gradient J^T e
	//	*                                 2 - stopped by small Dp
	//	*                                 3 - stopped by itmax
	//	*                                 4 - singular matrix. Restart from current p with increased \mu
	//	*                                 5 - no further error reduction is possible. Restart with increased mu
	//	*                                 6 - stopped by small ||e||_2
	//	*                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values; a user error
	//	* info[7]= # function evaluations
	//	* info[8]= # Jacobian evaluations
	//	* info[9]= # linear systems solved, i.e. # attempts for reducing error
	//	*/
	//	double *work,      /* I: working memory, allocated internally if NULL. If !=NULL, it is assumed to point to
	//					   * a memory chunk at least LM_DIF_WORKSZ(m, n)*sizeof(double) bytes long
	//					   */
	//	double *covar,     /* O: Covariance matrix corresponding to LS solution; Assumed to point to a mxm matrix.
	//					   * Set to NULL if not needed.
	//					   */
	//	void *adata)       /* I: pointer to possibly needed additional data, passed uninterpreted to func.
	//					   * Set to NULL if not needed
	//					   */
}

