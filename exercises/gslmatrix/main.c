#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// solve system Ax=b using Householder solver


void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n", s);
	for(int i=0; i<A->size1;i++){
		for(int j=0; j<A->size2; j++){
		printf("%10g",gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}



void vector_print(char s[], gsl_vector* u){
	printf("%s\n",s);
	for(int i=0; i<u->size; i++){
        	printf("%10g",gsl_vector_get(u,i));
        printf("\n");
	}
}


int main(){
int n = 3; //3x3 matrix
double a[] = {6.13, -2.90, 5.86, 8.08, -6.31, -3.89, -4.39, 1.00, 0.19}; //A
double b[] = {6.23, 5,27, 2.29}; //b

	gsl_matrix* A_copy = gsl_matrix_alloc(n,n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	gsl_matrix_view A = gsl_matrix_view_array(a,n,n);
	gsl_vector_view B = gsl_vector_view_array(b,n);
	gsl_matrix_memcpy(A_copy, &A.matrix);

	gsl_linalg_HH_solve(A_copy, &B.vector, x);
	gsl_blas_dgemv(CblasNoTrans, 1, &A.matrix, x, 0, y);

	vector_print("Righthand side b = ", &B.vector);
	printf("\n");
	vector_print("Check: A*x should be equal to b...:", y);
	printf("\n");

	gsl_matrix_free(A_copy);
	gsl_vector_free(x);
return 0;
}




