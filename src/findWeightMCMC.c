#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "concordanceIndex.h"

void findWeightMCMC(const double* x, const double* surv, double *w, int *mIn, int *nIn, double *dIn, int *maxIterIn){
	int i,j, count = 0, maxIter = *maxIterIn, m = *mIn, n = *nIn;
	double* pred = (double*) calloc(n, sizeof(double));
	double* bestPred = (double*) calloc(n, sizeof(double));
	double ccdi_cur, ccdi_max, delta = *dIn;


	Rprintf("Starting weights: \n");
	for(i = 0; i < m; i++) Rprintf("\t%f", w[i]);
	Rprintf("\n");

	for(j = 0; j < n; j++){
		bestPred[j] = 0;
		for(i = 0; i < m; i++){
			bestPred[j] += x[i + j*m] * w[i];
		}
	}

	
	Rprintf("Starting CCDI: \n");
	ccdi_max = concordance_index(bestPred, surv, n);
	Rprintf("\t%f", ccdi_max);
	Rprintf("\n");

	memcpy(pred, bestPred, n*sizeof(double));
	ccdi_cur = ccdi_max;

	//while(count < maxIter){
		

	//}

	free(bestPred);
	free(pred);
	Rprintf("Leave C function.\n");
}

