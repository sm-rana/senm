// test the class Varima


#include "Varima.h"
#include <stdio.h>

int main() {

	FILE* file = fopen("testVarima.dat", "w+");

	printf("Test 1: non-seasonal arima\n");
	Varima* varima = new Varima(3, 2, 0);

	double phi_in[] = {1.28, -0.4, 
					 0.7, 0,
					 1, -0.5};
	double cov[] = { 1, 0.3,  -0.2, 
					 0.3,  1,  0.5,
					 -0.2,  0.5, 1};
	varima->setPara(phi_in, NULL, cov);
	varima->initStateG(new Rvgs(1234), (double**)NULL, 0);
	
	double z_out[3];
	for (int istep=0; istep<100000; ++istep) {
		varima->generate(z_out);
		fprintf(file, "%f\t%f\t%f\n", z_out[0], z_out[1], z_out[2]);
	}
	fclose(file);

	////////////
	
	file = fopen("testVarima2.dat", "w+");

	printf("Test 2: double-seasonal arima s1=24, d=1\n");
	Varima* varima_diurnal = new Varima(varima, 1, 0, 24, 1);

	double phi_in2[] = {-0.4, 
					 0.7,
					 -0.12};
	varima_diurnal->setPara(phi_in2, NULL, NULL);
			
	for (int istep=0; istep<100000; ++istep) {
		varima_diurnal->generate(z_out);
		fprintf(file, "%f\t%f\t%f\n", z_out[0], z_out[1], z_out[2]);
	}

	fclose(file);
}