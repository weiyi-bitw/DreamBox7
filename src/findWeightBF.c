#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <time.h>
#include "concordanceIndex.h"

void findWeightBF(const double* x, const double* surv, double *w, int *mIn, int *nIn, double *dIn, int *maxIterIn, int *verbose){
	int i,j,k, count = 0, maxIter = *maxIterIn, m = *mIn, n = *nIn, streak = 0;
	double* dW = (double*) calloc(m-1, sizeof(double));
	double* pred = (double*) calloc(n, sizeof(double));
	double* bestPred = (double*) calloc(n, sizeof(double));
	double ccdi_cur, ccdi_max, delta_max = *dIn, ratio;

	srand(time(NULL));

	/*
	Rprintf("Starting weights: \n");
	for(i = 0; i < m; i++) Rprintf("\t%f", w[i]);
	Rprintf("\n");
	*/
	for(j = 0; j < n; j++){
		bestPred[j] = 0;
		for(i = 0; i < m; i++){
			bestPred[j] += x[i + j*m] * w[i];
		}
	}

	
	ccdi_max = concordance_index(bestPred, surv, n);
	Rprintf("Starting CCDI: \t%f\n", ccdi_max);

	memcpy(pred, bestPred, n*sizeof(double));
	//Rprintf("%d\t%d\n", count, maxIter);
	while(count < maxIter){
		for(k = 1; k < m; k++){
			//ratio = (double) (rand() % 100) / 100;
			ratio = (double) (rand() % 100 - 50) / 100;
			//dice = rand() % 3;
			if(rand() & 1) {
				//dW[k-1] = delta_max * ratio;
				dW[k-1] = ratio * w[k];
				for(j = 0; j < n; j++) pred[j] += dW[k-1] * x[(k) + j*m];
			//}else if(dice == 1){
				//dW[k-1] = -delta_max * ratio;
			//	dW[k-1] = ratio * w[k];
			//	for(j = 0; j < n; j++) pred[j] += dW[k-1] * x[(k) + j*m];
			}else{
				dW[k-1] = 0;
			}
			
		}
		ccdi_cur = concordance_index(pred, surv, n);
		//Rprintf("%f\n", ccdi_cur);
		if(ccdi_cur > ccdi_max){
			ccdi_max = ccdi_cur;
			for(k = 1; k < m; k++) w[k] += dW[k-1];
			
			if(*verbose == 1){
				Rprintf("New best CCDI: \t%f\t Streak: \t%d \n", ccdi_max, count);
				Rprintf("Weights: \n");
				for(i = 0; i < m; i++) Rprintf("\t%f", w[i]);
				Rprintf("\n");
			}

			if(count > streak) streak = count;
			count = 0;
			memcpy(bestPred,pred, n*sizeof(double));
		}else{
			memcpy(pred, bestPred, n*sizeof(double));
		}

		//if(count % 1000 == 0) Rprintf("%d\n", count);
		count++;
	}



	free(dW);
	free(bestPred);
	free(pred);
	
	Rprintf("================\n");
	Rprintf("Best CCDI found: \t%f \n", ccdi_max);
	Rprintf("Weights: \n");
	for(i = 0; i < m; i++) Rprintf("\t%f", w[i]);
	Rprintf("\n");
	Rprintf("\nLongest streak: \t%d \n", streak);
	Rprintf("================\n");
	//Rprintf("Leave C function.\n");
}

