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

/* Test classes Varima and Population. Multiple working threads are invoked 
to generate several 1-D vector series. Statistics are reported */

#include <stdio.h>
#include <Windows.h>
#include <process.h>
#include "Population.h"
#include "Varima.h"

//chain size
#define TEST_VARIMA_M  100 * 6  

struct Arg {
	unsigned iw;
	unsigned dim1;
	Varima* var;
	Population* pop;
};

unsigned WINAPI workVarimaGeneration(void* arg) {
	Arg* parg =  (Arg*) arg; //worker/chain id
	//double* z = new double[parg->dim1];
	double* z = (double*)calloc(parg->dim1, sizeof(double));

	parg->var->threadN = parg->iw;
	unsigned im, id;
	for (im = 0; im < TEST_VARIMA_M; ++im) {
		parg->var->generate(z);
		for (id = 0; id < parg->dim1; ++id) {
			POPH(parg->pop, parg->iw, im, id) = z[id];
		}
		//advance chain pointer
		parg->pop->p[parg->iw] ++;
        Sleep(30);
	}

	free(z);
    //delete[] z;
    return 0;
}




int main() {

	ewi(TEXT("Test 1: non-seasonal arima\n"));

	double phi_in[] = {
		1.28,	-0.4	, 
		0.7	,	0		,
		1	,	-0.5	};
	double cov[] = { 
		1,	0.3,	-0.2, 
		0.3,  1,	0.5,
		-0.2, 0.5,	1};

	Population* pop1;
	if (Pop_new(TEST_VARIMA_M, 3, &pop1)!=POP_OK) {return 1;}

    HANDLE worker[N_WORKERS];

	unsigned iw;
	for (iw=0; iw<N_WORKERS; ++iw) {
		Varima* var1 = new Varima(3, 2, 0);
		var1->setPara(phi_in, NULL, cov);
		var1->initStateG(new Rvgs(1234+iw), (double**)NULL, 0);

		Arg *parg = (Arg*)calloc(1,sizeof(Arg));
		parg->iw = iw;
		parg->pop = pop1;
		parg->var = var1;
		parg->dim1 = 3;

        worker[iw] = (void*)_beginthreadex(NULL, 0, 
			workVarimaGeneration, (void*)parg, 0, NULL);
		if (worker[iw] == NULL) {
			ewi(TEXT("Can't create thread.\n")); return 2;}
	}

    DWORD reason; //the reason why main thread is awaken
    do {
		reason = WaitForMultipleObjects(
			N_WORKERS, worker, TRUE, 4*1000); // wait 4 seconds
        ewi(TEXT("Calc stats...\n"));
        Pop_calc(pop1);
        Pop_report(pop1);
        fflush(stdout);
	} while (reason == WAIT_TIMEOUT);

    Pop_del(&pop1);

	//TODO: double-dimensional (d2 > 1) test 


    return 0;

}