#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_blas.h>


void matrixprint(gsl_matrix* A,FILE* file){
	int n=A->size1;
		int m=A->size2;

		for(int i=0;i<n;i++){
		        for(int j=0;j<m;j++){
		                fprintf(file,"%10.2g",gsl_matrix_get(A,i,j));
		        }
		        fprintf(file,"\n");}
}




void cholesky_decomp(gsl_matrix*A,gsl_matrix*L,int n){
//Cholesky–Crout algorithm.
//double eps = 0.01;
for(int j = 0; j < A->size1; j++){
	for(int i = j; i<A->size1; i++){
	if (i==j){
		double sum = 0;
		for (int k = 0; k < j; k++){
//			printf("k:%i\n",k);
			double Ljk = gsl_matrix_get(L,j,k);
			sum += pow(Ljk,2);
//			printf("jjsum:%g\n",sum);
		}
		double Ajj = gsl_matrix_get(A,j,j);
		double Ljj = sqrt(Ajj-sum);
		gsl_matrix_set(L,j,j,Ljj);
		}
	else{
//		printf("i j: %i %i\n",i,j);
		double sum = 0;
			for(int k = 0; k<j;k++){
//			printf("k:%i\n",k);
			double Ljk = gsl_matrix_get(L,j,k);
			double Lik = gsl_matrix_get(L,i,k);
//			printf("Ljk Lik: %g %g\n",Ljk,Lik);
			sum += (Lik*Ljk);
			}
//		printf("sum:%g\n",sum);
		double Ljj = gsl_matrix_get(L,j,j);
//		printf("Ljj:%g\n",Ljj);
		double Aij = gsl_matrix_get(A,i,j);
//		printf("Aij:%g\n",Aij);
		int Lij = ((Aij-sum)/Ljj);
//		printf("Lij:%g\n\n",Lij);
		gsl_matrix_set(L,i,j,Lij);
	}}
}
}




void vectorprint(gsl_vector* v, FILE* file){
	for(int i=0;i<v->size;i++)fprintf(file, "%10g ",gsl_vector_get(v,i));
	fprintf(file,"\n");
}

void backsub(gsl_matrix*U,gsl_vector*y,gsl_vector*c){
	int n = U->size1;
	for (int i =0; i<n; i++){
		double s = 0;
		for(int j = i+1; j<n ;j++){
			s += gsl_matrix_get(U,i,j)*gsl_vector_get(y,j);
		}
		gsl_vector_set(y,i,(gsl_vector_get(c,i)-s)/gsl_matrix_get(U,i,i));
	}


}
void forsub(gsl_matrix*L,gsl_vector*y,gsl_vector*c){
	int n = L->size1;
	for (int i =0; i<n; i++){
		double s = 0;
		for(int k = 1; k<i ;k++){
			s += gsl_matrix_get(L,i,k)*gsl_vector_get(y,k);
		}
		gsl_vector_set(y,i,(gsl_vector_get(c,i)-s)/gsl_matrix_get(L,i,i));
	}


}


void cholesky_solve(gsl_matrix*A, gsl_matrix*L, gsl_vector*b, gsl_vector*x,gsl_vector*y){
//	gsl_vector* y = gsl_vector_alloc(A->size1);

	gsl_blas_dgemv(CblasTrans,1.,L,x,0.,y);
	forsub(L,y,b);

	gsl_matrix*Lt = gsl_matrix_alloc(A->size1,A->size1);
	gsl_matrix_transpose_memcpy(Lt,L);
	backsub(Lt,x,y);



}







