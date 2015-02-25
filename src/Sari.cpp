#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "Sari.h"

Sari_Err Sari_new(unsigned n_in,  ///< dimension
				  unsigned ns_in, ///< seasonality
				  unsigned* H_in, ///< H matrix and size
				  unsigned H_len,
				  Sari** model_out) {

	if (H_len != 3 * (ns_in+1)) {
		return SARI_H_MATRIX_ERROR;
	}

	Sari* p = (Sari*)calloc(1, sizeof(Sari));
	if (p == NULL) return SARI_ALLOC_ERROR;

	p->H = (unsigned*)calloc(H_len, sizeof(unsigned));
	if (p->H == NULL) return SARI_ALLOC_ERROR;

	// copy n, ns, H
	p->n = n_in; p->ns = ns_in;
	memcpy(p->H, H_in, H_len*sizeof(unsigned));

	// compute look-back length
	unsigned lb = 0;
	unsigned phi_len = 0;
	for (unsigned i_ns = 0; i_ns <= p->ns; ++i_ns)  {
		// s * (p + d)
		lb += SARI_H(p, 0, i_ns)*(SARI_H(p, 1, i_ns)+SARI_H(p, 2, i_ns));
		phi_len += SARI_H(p, 2, i_ns);
	}
	p->lb = lb;
	p->phi_len = phi_len;

	// allocate mem for cov, phi, lcov
	p->cov = (double*)calloc(n_in * n_in, sizeof(double));
	p->lcov = (double*)calloc(n_in * n_in, sizeof(double));
	p->phi = (double*)calloc(phi_len, sizeof(double));
	p->poly = (double*)calloc(lb + 1, sizeof(double));
	if (p->cov == NULL || p->lcov == NULL || p->phi == NULL)
		return SARI_ALLOC_ERROR;

	*model_out = p;
	return SARI_OK;
}

void Sari_del(Sari** sari) {
	free((*sari)->H);
	free((*sari)->phi);
	free((*sari)->cov);
	free((*sari)->lcov);

	free(*sari);
}

///< multiply to a polynomial
static void poly_mlpl(double* poly_a, ///< a polynomial
					  unsigned* len_a,  ///< length of poly_a
					  double* poly_b, ///< another polynomial, 
									///< the first "1" is omitted
					  unsigned len_b, ///<length of poly_b
					  unsigned s ///< step of poly_b
					  ) {
	unsigned tp_len = *len_a + s * len_b;
	double* tp = (double*)calloc(tp_len, sizeof(double));

	memcpy(tp, poly_a, *len_a * sizeof(double));

	for (unsigned ib = 0; ib < len_b; ++ib) {
		for (unsigned ia = 0; ia < *len_a; ++ia) {
			tp[ia + (ib+1)*s] += poly_a[ia]*poly_b[ib];
		}
	}

	//update length of poly_a
	*len_a = tp_len;
	memcpy(poly_a, tp, tp_len * sizeof(double));
	free(tp);
}

Sari_Err Sari_setPhi(Sari* sari, double* phi_in, unsigned phi_len) {
	if (phi_len != sari->phi_len) return SARI_PHI_ARRAY_ERROR;
	memcpy(sari->phi, phi_in, sizeof(double)*phi_len);

	// compute poly[]: polynomial multiplication
	double diff_poly = -1.0;
	sari->poly[0] = 1;
	unsigned poly_len = 1;

	unsigned offset = 0;
	for (unsigned ins = 0 ; ins<=sari->ns; ++ins) {
		unsigned s = SARI_H(sari, 0, ins);
		unsigned p = SARI_H(sari, 2, ins);
		unsigned d = SARI_H(sari, 1, ins);
		poly_mlpl(sari->poly, &poly_len, &sari->phi[offset], p, s);
		for (unsigned id = 0; id < d; ++id) {
			poly_mlpl(sari->poly, &poly_len, &diff_poly, 1, s); 
		}
		offset += p;
	}

	// compute pos[] and para[]
	unsigned npos=0;
	unsigned ipoly=0; 
	for (; ipoly<=sari->lb; ++ipoly) {
		if (sari->poly[ipoly]!=0) {
			++npos;
		}
	}
	sari->npos = npos;

	sari->pos = (unsigned*)calloc(npos, sizeof(unsigned));
	sari->para = (double*)calloc(npos, sizeof(double));
	if (sari->pos == NULL || sari->para == NULL) return SARI_ALLOC_ERROR;

	unsigned ipos = 0;
	double ph = 0.0;
	for (ipoly=0; ipoly<=sari->lb; ++ipoly) {
		ph = sari->poly[ipoly]; // an AR parameter
		if (ph!=0) {
			sari->pos[ipos] = ipoly;
			sari->para[ipos] = ph;
			++ipos;
		}
	}

	return SARI_OK;
}

static void print_mat(unsigned* p_mat, ///<matrix to be printed, column-major
					  unsigned n_row, 
					  unsigned n_col, 
					  char* name) {
	if (p_mat == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s matrix (positive integers): (%u x %u)\n", 
		name, n_row, n_col);
	for (unsigned ir = 0; ir < n_row; ++ir) {
		for (unsigned ic = 0; ic < n_col; ++ic) {
			printf("%8u ", p_mat[ir + n_row * ic]);
		}
		printf("\n");
	}
}

static void print_mat(double* p_mat, ///<matrix to be printed, column-major
					  unsigned n_row, 
					  unsigned n_col, 
					  char* name) {
	if (p_mat == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s matrix (double-precision floats): (%u x %u)\n", 
		name, n_row, n_col);
	for (unsigned ir = 0; ir < n_row; ++ir) {
		for (unsigned ic = 0; ic < n_col; ++ic) {
			printf("%8.5f ", p_mat[ir + n_row * ic]);
		}
		printf("\n");
	}
}
static void print_vec(unsigned* p_vec, ///<vector to be printed
					  unsigned n, ///< number of elements 
					  char* name) {
	if (p_vec == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s vector (positive integers): (n = %u)\n", 
		name, n);
	for (unsigned i = 0; i < n; ++i) {
		printf("%8u ", p_vec[i]);
	}
	printf("\n");
}
static void print_vec(double* p_vec, ///<vector to be printed
					  unsigned n, ///< number of elements 
					  char* name) {
	if (p_vec == NULL) {
		printf("%s not allocated.\n", name);
		return;
	}
	printf("%s vector (double-precision floats): (n = %u)\n", 
		name, n);
	for (unsigned i = 0; i < n; ++i) {
		printf("%8.5f ", p_vec[i]);
	}
	printf("\n");
}



void Sari_dump(Sari* sari) {
	if (sari == NULL) {
		printf("Sari object not created.\n");
		return;
	}
	printf("Sari object %p, %u dimensions, %u layers of seasonality.\n",
		sari, sari->n, sari->ns);

	print_mat(sari->H, 3, sari->ns + 1, "H");
	print_vec(sari->phi, sari->phi_len, "phi");
	print_vec(sari->poly, sari->lb + 1, "poly");
	print_vec(sari->pos, sari->npos, "pos");
	print_vec(sari->para, sari->npos, "para");
	print_mat(sari->cov, sari->n, sari->n, "cov");
	print_mat(sari->cov, sari->n, sari->n, "lcov");

}








