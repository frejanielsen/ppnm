#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcs.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

int main(){
	int n = 3;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* L = gsl_matrix_alloc(n,n);
	gsl_matrix* A_c = gsl_matrix_alloc(n,n);


	double i1[]= {25,15,-5}; double i2[]={15,18,0}; double i3[]={-5,0,11};
	for(int j = 0; j<n;j++){
		gsl_matrix_set(A,j,0,i1[j]);
		gsl_matrix_set(A,j,1,i2[j]);
		gsl_matrix_set(A,j,2,i3[j]);
		}
/*
	unsigned int rand_seed=2;
	for (int i = 0; i<n; i++){
		for (int j = i; j<n; j++){
			double rand_n = 10.0*rand_r(&rand_seed)/RAND_MAX;
			gsl_matrix_set(A,j,i,rand_n);
			gsl_matrix_set(A,i,j,rand_n);
		}}*/
	FILE* out = fopen("out.txt","w");
	fprintf(out,"Matrix A:\n");
	matrixprint(A,out);

	fprintf(out,"Where A is a real symmetric positive-definite matrix.\n");



	fprintf(out,"\nCholesky decomposition is done on A.\n");
	gsl_matrix_memcpy(A_c,A);

	cholesky_decomp(A,L,n);
	fprintf(out,"Matrix A is decomposed to L\nL:\n");
	matrixprint(L,out);
	gsl_matrix* Lt = gsl_matrix_alloc(n,n); // allocate menory for
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,L,L,0,Lt);
	fprintf(out,"\nThe Cholesky decomposition has been succesfull as LL^T=A\nLL^T=\n");
	matrixprint(Lt,out);




	///////
	gsl_vector*x = gsl_vector_alloc(n);
	gsl_vector*b = gsl_vector_alloc(n);
	gsl_vector*y = gsl_vector_alloc(n);
	gsl_vector*bs = gsl_vector_alloc(n);
	fprintf(out,"\n\n\n");
	unsigned int rand_seed=2;
	for (int i = 0; i<n; i++){
	gsl_vector_set(b,i,1.0*rand_r(&rand_seed)/RAND_MAX);
	}
	//gsl_blas_dgemv(CblasTrans,1.,L,x,0.,y);


	fprintf(out,"Solving Ax=b using cholesky decomposition on matrix A\n");
	fprintf(out,"b=\n");
	vectorprint(b,out);
	cholesky_solve(A_c,L,b,x,y);

	fprintf(out,"\nL^T*x = y =\n");
	vectorprint(y,out);
	fprintf(out,"\nLy=b=\n");
	gsl_blas_dgemv(CblasNoTrans,1.,L,y,0,bs);
	vectorprint(bs,out);
return 0;
}
