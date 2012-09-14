#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

double concordance_index(const double* , const double*, R_len_t);
SEXP ccdiR2C(SEXP, SEXP);
void concordance_index_all(const double*, const double*, const int*, const int*, double*);
