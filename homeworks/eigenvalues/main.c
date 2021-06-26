#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <time.h>

#include "evfs.h"


int main(){
int n = 6; //size of matricies

gsl_matrix* A = gsl_matrix_alloc(n,n);
gsl_matrix* V = gsl_matrix_alloc(n,n);
gsl_matrix* A_copy = gsl_matrix_alloc(n,n);
gsl_matrix* Prod = gsl_matrix_alloc(n,n);
gsl_matrix* Prod_ = gsl_matrix_alloc(n,n);

// fill A with random numbers. V is identity matrix
double r=0;
unsigned int rand_seed=1;
for (int i=0; i<n;i++){
	for (int j=i;j<n;j++){
		r=1.0*rand_r(&rand_seed)/RAND_MAX;
		gsl_matrix_set(A,i,j,r);
		gsl_matrix_set(A,j,i,r);
		gsl_matrix_set(V,i,j,0);
	}
	gsl_matrix_set(V,i,i,1);
}

// copy A
gsl_matrix_memcpy(A_copy,A);

// print matricies
printf("Before diagonalisation:\n A=\n");
matrixprint(A);

printf("V=\n");
matrixprint(V);

jacobi_diag(A,V);
printf("After diagonlisation:\n");
printf("D=\n");
matrixprint(A);

printf("V=\n");
matrixprint(V);

printf("V^T A V=\n");
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_copy,V,0.0,Prod);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,Prod,0.0,Prod_);
matrixprint(Prod_);


printf("V D V^T=\n");
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,V,0.0,Prod);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,Prod,0.0,Prod_);
matrixprint(Prod_);

printf("V^T V=\n");
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,V,0.0,Prod);
matrixprint(Prod);

////////// Quantum particle in a box /////////////

int N = 20;
double s=1.0/(N+1);
gsl_matrix* H = gsl_matrix_alloc(N,N);

for (int i=0;i<N-1;i++){
	gsl_matrix_set(H,i,i,-2);
	gsl_matrix_set(H,i,i+1,1);
	gsl_matrix_set(H,i+1,i,1);
	}

gsl_matrix_set(H,N-1,N-1,-2);
gsl_matrix_scale(H,-1/s/s);

gsl_matrix* Psi = gsl_matrix_alloc(N,N); // Psi -> J
for (int i=0;i<N;i++){
	gsl_matrix_set(Psi,i,i,1);
}

// diagonalize H
printf("Before diagonalization:\n H=\n");
matrixprint(H);

printf("States=\n");
matrixprint(Psi);

jacobi_diag(H,Psi);

printf("After diagonalization:\n H=\n");
matrixprint(H);
printf("States=\n");
matrixprint(Psi);

// check energies

printf("Energies:\n");
printf("Nummerical   Exact\n");

for (int k=0; k<N;k++){
	double exact_E = M_PI*M_PI*(k+1)*(k+1);
	double calc_E = gsl_matrix_get(H,k,k);
	printf("%i %g %g\n",k ,calc_E ,exact_E);
	}
FILE* States_file = fopen("out_states.txt","w");
for (int i=0; i<3; i++){
	fprintf(States_file,"%10g %10g\n",0.0,0.0);
	for (int j=0;j<N;j++){
		fprintf(States_file,"%10g %10g\n",1.0*(j+1)/(N+1),gsl_matrix_get(Psi,j,i));
	}
	fprintf(States_file,"%10g %10g",1.0,0.0);
	fprintf(States_file,"\n\n\n\n\n");
}






return 0;
}
