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

	printf("Q^t*Q=\n");
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,prod);
	matrixprint(prod);

	printf("we conclude the QR-decomposition works for tall matricies \n");

//	printf("test \n");

///////// fitting ///////

	int N_f = 2;
	gsl_vector* c = gsl_vector_alloc(N_f);

	FILE* data_file2 = fopen("out_data.txt","w");

	for (int i=0; i<N; i++){
	        for (int j=0; j<2; j++){
	                fprintf(data_file2,"%10g", gsl_matrix_get(data,i,j));
	        }
		fprintf(data_file2," %10g \n",gsl_matrix_get(data,i,1)/20);
	}

	       	for (int j=0; j<N; j++){
	           	gsl_matrix_set(data,j,1,log(gsl_matrix_get(data,j,1)));
			gsl_matrix_set(data,j,2,1./20);
		}

	printf("\n");
	//printf("The data we will be fitting to:\n");
	matrixprint(data); //  print data matrix


	gsl_matrix* cov = gsl_matrix_alloc(N_f,N_f); //covariance
	fit(data, &func, N_f, c, cov);
	printf("covariance matrix =\n");
	matrixprint(cov);

	printf("C=\n");

	       for (int j=0; j<N_f; j++){
	                printf("%10g   ",gsl_vector_get(c,j));
	        }

	printf("\nDelta C=\n");

	       for (int j=0; j<N_f; j++){
	                printf("%10g   ",sqrt(gsl_matrix_get(cov,j,j)));
	        }


	printf("\n");



	double mht = log(1./2)/gsl_vector_get(c,1);
	double mhtM = log(1./2)/(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(cov,1,1)));
	double mhtm = log(1./2)/(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(cov,1,1)));

	printf("\n Halftime calculated:\n %10g Min: %10g max: %10g",mht,mhtM,mhtm);
	printf("\n Online: 3.63\n the two dont agree."); // However it should benoted that the other isotopes of Radium has markedly higher half times (14 days and 1600 years respectively) and so contamination could lead to the higher half-time observed.");

	FILE* fit_stream = fopen("out_fit.txt","w");

	       for (int j=1; j<=15; j++){
	                fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-gsl_vector_get(c,1)*j));
	        }
	 	fprintf(fit_stream,"\n\n\n\n\n\n\n\n");
	      for (int j=1; j<=15; j++){
	          fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)+sqrt(gsl_matrix_get(cov,0,0))-(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(cov,1,1)))*j));
	      }
	        fprintf(fit_stream,"\n\n\n\n\n\n\n\n");

		for (int j=1; j<=15; j++){
	          fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-sqrt(gsl_matrix_get(cov,0,0))-(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(cov,1,1)))*j));
	      }
	        printf("\n\n\n\n\n\n\n\n");



gsl_matrix_free(cov);
gsl_vector_free(c);
gsl_matrix_free(prod);
gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(data);

}
