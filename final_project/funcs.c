#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_blas.h>

void cholesky_decomp(gsl_matrix*A,gsl_matrix*L){
	int size_c = A->size2; // number of cols
	for (int i=0; i<size_c;i++){
		for(int j=0; j<=i;j++){
		double sum = 0;
		for (int k=0; k<j; j++){
			double L_ik = gsl_matrix_get(L,i,k);
			double L_jk = gsl_matrix_get(L,j,k);
			sum += (L_ik + L_jk);
		}
		if(i==j){
			double Lij = sqrt(gsl_matrix_get(A,i,i)- sum);
			gsl_matrix_set(L,i,j,Lij);
		}
		else{
			double Lij = (1.0/gsl_matrix_get(L,j,j)*(gsl_matrix_get(A,i,j)-sum));
			gsl_matrix_set(L,i,j,Lij);
		}
	}
}}


void matrixprint(gsl_matrix* A,FILE* file){
	int n=A->size1;
		int m=A->size2;

		for(int i=0;i<n;i++){
		        for(int j=0;j<m;j++){
		                fprintf(file,"%10.2g",gsl_matrix_get(A,i,j));
		        }
		        fprintf(file,"\n");}
}

/*
double randomNumber( unsigned int *seed ){     //generate random numbers to test program.
	double maxRand = (double)RAND_MAX;
	double randNum = (double)rand_r( seed );
	return randNum/maxRand;
}


void vectorprint(gsl_vector* v, FILE* file){
	for(int i=0;i<v->size;i++)fprintf(file, "%10g ",gsl_vector_get(v,i));
}

void cholesky_solve(gsl_matrix*A, gsl_matrix*L, gsl_vector*b, gsl_vector*x){





}






*/

