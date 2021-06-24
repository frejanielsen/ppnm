#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "ls_fs.h"

int main(){

	int n = 7;
	int m = 6;

	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_matrix* prod = gsl_matrix_alloc(m,m);


	for(int i=0; i<n; i++){ //Creates a random matrix of size n,m
	        for(int j=0; j<m; j++){
	            double A_ij = rand()/100;
	            gsl_matrix_set(A,i,j,A_ij);
	        }
	    }

	FILE* data_file = fopen("data.txt","r");
	int N = 9;  // number of data points
	int i = 0;

	gsl_matrix* data = gsl_matrix_alloc(N,3);
	double X, Y;
	while(2== fscanf(data_file, "%lg %lg", &X,&Y)){

	gsl_matrix_set(data,i,0,X);
	gsl_matrix_set(data,i,1,Y);
	i++;
	}

	printf("X-data: Y-data:\n");
	for(i=0;i<N;i++){
	printf("%10g %10g\n",gsl_matrix_get(data,i,0),gsl_matrix_get(data,i,1));}


	printf("A=\n");
	matrixprint(A);

	GS_decomp(A,R);
	printf("Q=\n");

	matrixprint(A);

	printf("R=\n");
	matrixprint(R);




	printf("test \n");


}
