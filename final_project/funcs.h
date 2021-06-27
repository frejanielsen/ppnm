#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

void cholesky_decomp(gsl_matrix*A,gsl_matrix* L);
void matrixprint(gsl_matrix* A,FILE* file);
