#include "splineMI.h"

void quickSortR(float const * arr, int *idx, int left, int right) {
      int i = left, j = right;
      int tmp;
      float pivot = arr[ idx[(left + right)/2] ];

      /* partition */
      while (i <= j) {
            while (arr[idx[i]] < pivot)
                  i++;
            while (arr[idx[j]] > pivot)
                  j--;
            if (i <= j) {
                  tmp = idx[i];
                  idx[i] = idx[j];
                  idx[j] = tmp;
                  i++;
                  j--;
            }
      };

      if (left < j)
            quickSortR(arr, idx, left, j);
      if (i < right)
            quickSortR(arr, idx, i, right);
}

double entropy1s(const double *weights, const bool *eff, int numSamples, int numBins) {
  int curBin, curSample, k = numSamples;
  double H = 0, h;
  for(curSample = 0; curSample < numSamples; curSample++){
	if(!eff[curSample]) k--;
  }

  for (curBin = 0; curBin < numBins; curBin++) {
    h = 0;
    for(curSample = 0; curSample < numSamples; curSample++){
	if(eff[curSample]) h += weights[curBin * numSamples + curSample];
    }
    h /= k;
    if (h > 0) {
      H -= h * log2d(h);
    }
  }
  return H;
}

double entropy2sDiffBins(const double *wx, const double *wy, const bool *eff, int numSamples, int binx, int biny){
	int curBinX, curBinY, curSample, k = numSamples;
	double H = 0, h;
	for(curSample = 0; curSample < numSamples; curSample++){
		if(!eff[curSample]) k--;
	}
	for(curBinX = 0; curBinX < binx; curBinX++){
	  for(curBinY = 0; curBinY < biny; curBinY++){
		h = 0;
		for(curSample = 0; curSample < numSamples; curSample++){
			if(eff[curSample]) h += wx[curBinX*numSamples + curSample] * wy[curBinY * numSamples + curSample];
		}
		h /= k;
		if(h > 0){
			H -= h * log2d(h);
		}
	  }
	}
	return H;
}


/*
void order(float const *x, int n, int *o){
  int i;
  for(i = 0; i < n; i++){
    o[i] = i;
  }
  quickSortR(x, o, 0, n-1);
}
*/

void msmi(double const *x, double const *surv, int *nIn, int *binx, int *sox, double *tstart, double *tend, double *misum, double *mimax, double *tmax){
	int i,ti, n = *nIn;
	int *o = (int*) calloc(n, sizeof(int));
	double *ux = (double*) calloc(*binx + *sox, sizeof(double));
	double *wx = (double*) calloc(*binx * n, sizeof(double));
	bool *effective = (bool*) calloc(n, sizeof(bool));
	double *wb = (double*) calloc(2 * n, sizeof(double));
	double e1x, e1b, curMI;

	*misum = 0; *mimax = 0; *tmax = 0;

	for(i = 0; i < n; i++) o[i] = i;
	quickSortR(surv, o, 0, n-1);

	knotVector(ux, *binx, *sox);
	findWeights(x, ux, wx, n, *sox, *binx, -1, -1);
	
	//first round find the first sample > tstart
	ti = 0;
	while(surv[o[ti]] < *tstart) ti++;
	for(i = 0; i < n; i++){
		if(i <= ti){
			wb[o[i]] = 1; wb[n + o[i]] = 0;
			effective[o[i]] = surv[n + o[i]] == 1? true : false;
		}else{
			wb[o[i]] = 0; wb[n + o[i]] = 1;
			effective[o[i]] = true;
		}
	}
	e1x = entropy1s(wx, effective, n, *binx)
	e1b = entropy1s(wb, effective, n, 2);
	curMI = e1x + e1b - entropy2sDiffBins(wx, wb, effective, n, *binx, 2);
	*mimax = curMI;
	*tmax = surv[o[ti]];
	*misum = curMI * (surv[o[ti]] - *tstart);
	ti++;
	while(surv[o[ti+1]] < *tend){
		wb[o[ti]] = 1; wb[n + o[i]] = 0;		
		effective[o[ti]] = surv[n + o[i]] == 1? true : false;
		e1x = entropy1s(wx, effective, n, *binx)
		e1b = entropy1s(wb, effective, n, 2);
		curMI = e1x + e1b - entropy2sDiffBins(wx, wb, effective, n, *binx, 2);
		if(curMI > *mimax){
			*mimax = curMI;
			*tmax = surv[o[ti]];
		}
		*misum = curMI * (surv[o[ti]] - surv[o[ti-1]]);

		ti++;
	}

	free(o);
	free(ux);
	free(wx);
	free(effective);
	free(wb);
}
