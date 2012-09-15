#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <time.h>
#include "concordanceIndex.h"

void findWeightMCMC(const double* x, const double* surv, double *w, int *mIn, int *nIn, int *kmaxIn, double *dIn, int *maxIterIn){
	int i,j,k, count = 0, maxIter = *maxIterIn, m = *mIn, n = *nIn;
	int kMax = *kmaxIn; // number of changed features
	double* dW = (double*) calloc(kMax, sizeof(double));
	int* dWi = (int*) calloc(kMax, sizeof(int));
	double* pred = (double*) calloc(n, sizeof(double));
	double* bestPred = (double*) calloc(n, sizeof(double));
	double ccdi_cur, ccdi_max, delta = *dIn, diff, ratio;

	srand(time(NULL));


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
	//Rprintf("%d\t%d\n", count, maxIter);
	while(count < maxIter){
		for(k = 0; k < kMax; k++){
			i = rand() % (m-1);
			dWi[k] = i+1;
			if(rand() & 1) {
				for(j = 0; j < n; j++) pred[j] += delta * x[(i+1) + j*m];
				dW[k] = delta;
			}else{
				for(j = 0; j < n; j++) pred[j] -= delta * x[(i+1) + j*m];
				dW[k] = -delta;
			}
		}
		ccdi_cur = concordance_index(pred, surv, n);
		diff = 100 * (ccdi_cur-0.5) / (ccdi_max-0.5);
		//Rprintf("%f\n", ccdi_cur);
		if(rand() % 100 <= diff){
			Rprintf("New best CCDI: \t%f \n", ccdi_cur);
			Rprintf("Weights: \n");
			for(i = 0; i < m; i++) Rprintf("\t%f", w[i]);
			Rprintf("\n");
			count = 0;
			ccdi_max = ccdi_cur;
			for(k = 0; k < kMax; k++) w[dWi[k]] += dW[k];
			memcpy(bestPred,pred, n*sizeof(double));
		}else{
			memcpy(pred, bestPred, n*sizeof(double));
		}


		count++;
	}


	free(dWi);
	free(dW);
	free(bestPred);
	free(pred);
	Rprintf("Leave C function.\n");
}

