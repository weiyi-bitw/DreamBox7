#include <stdio.h>

void concordance_index(double* predictions, double* observations, int* nIn, double* c)
{
	int count = 0, score = 0, i, j, n = *nIn;
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
	*c = (double) score/count;
}

double concordance_index_all(double* predictions, double* observations, int* mIn, int* nIn, double *score)
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

