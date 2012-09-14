#include "concordanceIndex.h"

double concordance_index(const double* predictions, const double* observations, R_len_t n)
{
	int count = 0, score = 0, i, j;
	double c;
	for(i = 0; i < n-1; i++){
	  for(j = i+1; j < n; j++){
		if( ((observations[n + i]==1) && (observations[n + j]==1)) ||
		((observations[n + i] == 1) && (observations[j] > observations[i])) || 
		((observations[n + j] == 1) && (observations[i] > observations[j])) )
		{
			count += 2;
			if(observations[i] < observations[j]){
				if(predictions[i] > predictions[j]) score += 2;
				else if(predictions[i] == predictions[j]) score++;
			}
			if(observations[i] == observations[j] && predictions[i] == predictions[j]) score += 2;
			if(observations[i] > observations[j]){
				if(predictions[i] < predictions[j]) score += 2;
				else if(predictions[i] == predictions[j]) score++;
			}


		}
	  }
	}
	c = (double) score/count;
	return c;
}

SEXP ccdiR2C(SEXP predIn, SEXP obsIn){
	SEXP out;
	double *predictions, *observations, *o;
	R_len_t n;

	PROTECT(predIn = AS_NUMERIC(predIn));
	PROTECT(obsIn = AS_NUMERIC(obsIn));
	PROTECT(out = NEW_NUMERIC(1));
	n = LENGTH(predIn);

	predictions = NUMERIC_POINTER(predIn);
	observations = NUMERIC_POINTER(obsIn);
	o = NUMERIC_POINTER(out);

	*o = concordance_index(predictions, observations, n);
	UNPROTECT(3);
	return out;

}

void concordance_index_all(const double* predictions, const double* observations, const int* mIn, const int* nIn, double *score)
{
        int count = 0, m = *mIn, n = *nIn, i, j, k;
        // initialize score array
        for(k = 0; k < m; k++) score[k] = 0;

        for(i = 0; i < n-1; i++){
          for(j = i+1; j < n; j++){
                if( ((observations[n + i]==1) && (observations[n + j]==1)) ||
                ((observations[n + i] == 1) && (observations[j] > observations[i])) ||
                ((observations[n + j] == 1) && (observations[i] > observations[j])) ){
                        count += 2;
                        if(observations[i] < observations[j]){
                                for(k = 0; k < m; k++){
                                        if(predictions[i * m + k] > predictions[j * m  + k]) score[k] += 2;
                                        else if(predictions[i * m + k] == predictions[j * m + k]) score[k]++;
                                }
                        }
                        if(observations[i] == observations[j]){
                                for(k = 0; k < m; k++){
                                        if(predictions[i*m + k] == predictions[j*m + k]) score[k] += 2;
                                }
                        }
                        if(observations[i] > observations[j]){
                                for(k = 0; k < m; k++){
                                        if(predictions[i*m + k] < predictions[j*m + k]) score[k] += 2;
                                        else if(predictions[i*m + k] == predictions[j*m + k]) score[k]++;
                                }
                        }



                }

          }
        }
        for(k = 0; k < m; k++) score[k] = score[k]/count;

}



/*
int main(){
	double p[7] = {1, 3, 2, 7, 8, 6, 4};
	double s[14] = {7, 6, 5, 4, 3, 2, 1, 1, 1, 0, 0, 1, 1, 1};

	double c = concordance_index(p, s, 7);

	printf("%f\n", c);
	return 0;
}
*/

