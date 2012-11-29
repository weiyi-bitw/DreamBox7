#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

double concordance_index(const double* , const double*, R_len_t);
SEXP ccdiR2C(SEXP, SEXP);
double concordance_index_w(const double* , const double*, const double*, R_len_t);
SEXP ccdiwR2C(SEXP, SEXP, SEXP);
void concordance_index_all(const double*, const double*, const int*, const int*, double*);
