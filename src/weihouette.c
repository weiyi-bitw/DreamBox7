#include <R.h>

void lowerTriIndex(int *x, int *y, int *m, int *diag){
	int i, j,c = 0;

	if(*diag == 1){
		for(i = 0; i < *m; i++){
			for(j = i; j < *m; j++){
				x[c] = i;
				y[c] = j;
				c++;
			}
		}
	}else{
		for(i = 0; i < *m; i++){
			for(j = (i+1); j < *m; j++){
				x[c] = i;
				y[c] = j;
				c++;
			}
		}
	}
}

void weihouette(double const *mis, int const *clust_member, int *k, int *m, int *lenpw, double *w, double *a, double *b){
	int i,j,h, numWithin = 0, diag, oi;
	double denom;
	//Rprintf("calloc...\n");
	int *x = (int*) calloc(*lenpw, sizeof(int));
	int *y = (int*) calloc(*lenpw, sizeof(int));
	double *crossmi = (double*) calloc(*k-1, sizeof(double));
	int *numCross = (int*) calloc(*k-1, sizeof(int));

	//Rprintf("diag[0] = 0...\n");
	diag = 0;
	
	//Rprintf("lowerTriIndex...\n");
	lowerTriIndex(x, y, m, &diag);
	
	for(i = 0; i < *k; i++){
		a[i] = 0;
		b[i] = 0;
		w[i] = 0;
		numWithin = 0;
		for(h = 0; h < *k-1; h++){
			numCross[h] = 0;
			crossmi[h] = 0;
		}
		Rprintf("k = %d ... ", i);
		for(j = 0; j < *lenpw; j++){
			if(clust_member[x[j]] == i && clust_member[y[j]] == i){
				numWithin ++;
				a[i] += mis[j];
			}else if(  clust_member[x[j]]==i || clust_member[y[j]] == i   ){
				oi = clust_member[x[j]]==i ? clust_member[y[j]] : clust_member[x[j]];
				if(oi > i) oi--;
				numCross[oi]++;
				crossmi[oi] += mis[j];
			}
		}
		for(h = 0; h < *k-1; h++){ 
			Rprintf("\t %f (%d)", crossmi[h], numCross[h]);
			crossmi[h] /= numCross[h];
			if(crossmi[h] > b[i]) b[i] = crossmi[h];
		}
		Rprintf("\n");
		a[i] /= numWithin;
		denom = a[i] < b[i] ? 1 - a[i] : 1-b[i];
		w[i] = (a[i] - b[i]) / denom;
	}
	free(x);
	free(y);
	free(crossmi);
	free(numCross);
}


