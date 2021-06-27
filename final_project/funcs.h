#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

void cholesky_decomp(gsl_matrix*A,gsl_matrix* L, int n);
void matrixprint(gsl_matrix* A,FILE* file);
void vectorprint(gsl_vector*v, FILE*file);
void backsub(gsl_matrix*U,gsl_vector*y,gsl_vector*c);
void forsub(gsl_matrix*L,gsl_vector*y,gsl_vector*c);
void cholesky_solve(gsl_matrix*A,gsl_matrix*L,gsl_vector*b,gsl_vector*x,gsl_vector*y);
