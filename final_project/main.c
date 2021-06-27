#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcs.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>


int main(){
	int n = 5;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* L = gsl_matrix_alloc(n,n);
	gsl_matrix* A_c = gsl_matrix_alloc(n,n);

	unsigned int rand_seed=2;
	for (int i = 0; i<n; i++){
		for (int j=0; j<n; j++){
			gsl_matrix_set(A,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		}}
	FILE* out = fopen("out.txt","w");
	fprintf(out,"A:\n");
	matrixprint(A,out);




/*
	int n = 4;
	double m1[n*n];
	srand(0);
	for (int i = 0; i < (3*n); i++) {
    		m1[i] = rand();
	}
	double *c1 = cholesky(m1, n);
	printf("m1:\n");
	print_matrix(c1, n);
	printf("\n");
	free(c1);
	n = 4;
	printf("m2:\n");
	double m2[16] = {18, 22,  54,  42,
	               22, 70,  86,  62,
	               54, 86, 174, 134,
	               42, 62, 134, 106};
	double *c2 = cholesky(m2, n);
	double* C2 = Cholesky(m2,n);
	print_matrix(c2, n);
	print_matrix(C2,n);
	free(c2);
	free(C2);
*/
return 0;
}
