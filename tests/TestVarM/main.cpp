#include "stdio.h"
#include "VarModel.h"

int main() {

/*Two test cases:
Case 1:
  n = 2;
  s = [ 1 24 ];
  p = [ 2, 1];
  phi1 = [ -.8  0.4;
           0.3  0.6];
  phi2 = [ 0.2  0; 
		   0.1  -0.2];
  phi3 = [ 0.9  -0.4;
		   -0.6  0.7];
  cov = [ 1  0;
          0  2];

Case 2:
  n = 2;
  s = [ 1 3 24];
  p = [ 3 1 1];
  phi1, phi2, phi3, same as above,
  phi4 = [0.4 0.3;
		  0.2 0.6];
  phi5 = [0.7 0; 
          0  0.7];
  cov = [ 1  0.5;
          0.5  1];

Case 3:
  n = 2;
  s = [ 1 ];
  p = [ 1 ];
  phi1 = [0.4 0.3;
		  0.2 0.6];
  cov = [ 1  0.5;
          0.5  1];
  panel data:
    y1 =	0.3		0.6		-0.5
	y2 =	0		-0.2	 0.5 

*/
	int s1[] = {1, 24};
	int p1[] = {2, 1};

	VarModel_Err err;
	VarModel* vm;
	err = VarModel_new(2, 1, s1, p1, &vm);
	if (err != VARMODEL_OK) {
		printf("Problem Creating VarModel 1.\n");
		return 1;
	}

	double phi1[] = {-0.8, 0.3, 0.4, 0.6, 
					 0.2, 0.1, 0, -0.2,
					 0.9, -0.6, -0.4, 0.7};
	err = VarModel_setPhi(vm, phi1, 12);
	if (err != VARMODEL_OK) {
		printf("Problem Setting VarModel 1 Phi parameters.\n");
		return 2;
	}

	double cov1[] = {1, 0, 0, 2};
	err = VarModel_setCov(vm, cov1, 4);
	if (err != VARMODEL_OK) {
		printf("Problem Setting VarModel 1 covariance matrix.\n");
		return 3;
	}

	VarModel_dump(vm);
	VarModel_del(&vm);

	//Case 2

	int s2[] = {1,3,24};
	int p2[] = {3,1,1};

	VarModel_Err err2;
	VarModel* vm2;
	err2 = VarModel_new(2, 2, s2, p2, &vm2);
	if (err2 != VARMODEL_OK) {
		printf("Problem Creating VarModel 2.\n");
		return 1;
	}

	double phi2[] = {-0.8, 0.3, 0.4, 0.6, 
					 0.2, 0.1, 0, -0.2,
					 0.9, -0.6, -0.4, 0.7,
					 0.4, 0.2, 0.3, 0.6,
					 0.7, 0, 0, 0.7};
	err2 = VarModel_setPhi(vm2, phi2, 20);
	if (err2 != VARMODEL_OK) {
		printf("Problem Setting VarModel 2 Phi parameters.\n");
		return 2;
	}

	double cov2[] = {1, 0.5, 0.5, 1};
	err2 = VarModel_setCov(vm2, cov2, 4);
	if (err2 != VARMODEL_OK) {
		printf("Problem Setting VarModel 2 covariance matrix.\n");
		return 3;
	}


	VarModel_dump(vm2);
	VarModel_del(&vm2);

	//Case 3

	int s3[] = {1};
	int p3[] = {1};

	VarModel_Err err3;
	VarModel* vm3;
	err3 = VarModel_new(2, 0, s3, p3, &vm3);
	if (err3 != VARMODEL_OK) {
		printf("Problem Creating VarModel 3.\n");
		return 1;
	}

	double phi3[] = { 0.4, 0.2, 0.3, 0.6 };
	err3 = VarModel_setPhi(vm3, phi3, 4);
	if (err3 != VARMODEL_OK) {
		printf("Problem Setting VarModel 3 Phi parameters.\n");
		return 2;
	}

	double cov3[] = {1, 0.5, 0.5, 1};
	err = VarModel_setCov(vm3, cov3, 4);
	if (err != VARMODEL_OK) {
		printf("Problem Setting VarModel 2 covariance matrix.\n");
		return 3;
	}

	//double panel[] = {0.3, 0, 0.6, -0.2};
	double panel[] = {0.3, 0, 0.6, -0.2, -0.5, 0.5};
	double llh;
	int nea;
	err = VarModel_logLikelihood(vm3, panel, 2, 3, &llh, &nea);
	if (err != VARMODEL_OK) {
		printf("Problem computing log-likelihood.\n");
		return 4;
	}
	printf("log-likelihood: %8.4f\n", llh);


	VarModel_dump(vm3);
	VarModel_del(&vm3);







}

