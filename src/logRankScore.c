#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void quickSortSurv(double const * arr, int const * event, int *idx, int left, int right) {
      int i = left, j = right;
      int tmp;
      double pivot = arr[ idx[(left + right)/2] ];
      int pivotEvent = event[ idx[(left + right)/2]];

      /* partition */
      while (i <= j) {
            while (arr[idx[i]] < pivot || ( (arr[idx[i]]==pivot) && (pivotEvent < event[idx[i]]) ))
                  i++;
            while (arr[idx[j]] > pivot || ( (arr[idx[j]]==pivot) && (pivotEvent > event[idx[j]]) ))
                  j--;
            if (i <= j) {
                  tmp = idx[i];
                  idx[i] = idx[j];
                  idx[j] = tmp;
                  i++;
                  j--;
            }
      };

	//printf("i:%d\tj:%d\tleft:%d\tright:%d\n", i, j, left, right);
      if (left < j)
            quickSortSurv(arr, event, idx, left, j);
      if (i < right)
            quickSortSurv(arr, event, idx, i, right);
}
/*
void orderSurv(double const *x, int const *event, int n, int *o){
  int i;
  for(i = 0; i < n; i++){
    o[i] = i;
  }
  quickSortSurv(x, event, o, 0, n-1);
}
*/
void logRankScore(double const *x, int const *evt, int *nIn, int* o, double *lrScore){
	int i, j, tieStart = 0, m = 0, n = *nIn;
	double preScore = 0, scoreCum = 0, thisScore;
	double *sortedScore = (double*) calloc(n, sizeof(double));
	
	for(i = 0; i < n; i++){
		o[i] = i;
		//REprintf("x:%f\to:%d\n", x[i], o[i]);
	}
	
	quickSortSurv(x, evt, o, 0, n-1);
	//for(i = 0; i < n; i++) REprintf("x:%f\to:%d\n", x[o[i]], o[i]);
	if(evt[o[0]]==1){
		preScore = (double)1/n;
		scoreCum = preScore-1;
		m = 1;
		tieStart = 0;
	}else{
		sortedScore[0] = 0;
	}
	for(i = 1; i < n; i++){
		if(evt[o[i]] == 1){
			preScore += (double)1 / (n-i);
			thisScore = preScore - 1;
			if( (x[o[i-1]] == x[o[i]]) && (evt[o[i-1]] == 1) ){
				scoreCum += thisScore;
				m++;
				if(i==n-1) for(j = tieStart; j < n; j++) sortedScore[j] = scoreCum/m;
			}else{
				for(j = tieStart; j < i; j++) sortedScore[j] = scoreCum/m;
				tieStart = i;
				m = 1;
				scoreCum = thisScore;
				if(i==n-1) sortedScore[i] = thisScore;
			}
			//printf("%f\n", preScore);
		}else{
			for(j = tieStart; j < i; j++) sortedScore[j] = scoreCum/m;
			tieStart = i+1;
			m = 0;
			scoreCum = 0;
			sortedScore[i] = preScore;
		}
	}
	
	for(i = 0; i < n; i++){
		lrScore[o[i]] = sortedScore[i];
	}

	free(sortedScore);
}

/*
int main(){
	double x[20] = {1, 3, 8, 9, 11, 13, 15, 20, 20, 3, 22, 25, 25, 31, 10, 33, 35, 40, 45, 30};
	int evt[20] = {1, 1, 1, 1,  1,  0,  1,  1,  0, 1,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0};
	double *lr = (double*) calloc(20, sizeof(double));
	int *o = (int*) calloc(20, sizeof(int));
	int i;
	logRankScore(x, evt, 20, o, lr);

	for(i = 0; i < 20; i++){
		printf("%f\t%d\t%f\n", x[o[i]], evt[o[i]], lr[o[i]]);
	}
	
	free(o);
	free(lr);
	return 0;
}
*/
