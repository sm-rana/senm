#include "Sari.h"

int main() {
	/* H matrix tested:
	 H1 = 
		1	24
		0	0
		1	1
	  phi0 = -0.5
	  phi1 = 0.3

	 H2 =
		1	24	168
		0	1	1
		1	3	1
	  phi0 = 0.3
	  phi1 = 0.5, -0.3, 0.1
	  phi2 = 0.4
		*/
	   
	unsigned H1[] = {1,0,1,24,0,1};
	double phi1[] = {-0.5,0.3};
	unsigned H2[] = {1,0,1,24,1,3,168,1,1};
	double phi2[] = {0.3,0.5,-0.3,0.1,0.4};

	Sari_Err err;
	Sari *sari;
	err = Sari_new(2, 1, H1, 3*2, &sari);
	err = Sari_setPhi(sari, phi1, 2);
	Sari_dump(sari);
	Sari_del(&sari);

	err = Sari_new(2, 2, H2, 3*3, &sari);
	err = Sari_setPhi(sari, phi2, 5);
	Sari_dump(sari);
	Sari_del(&sari);




}