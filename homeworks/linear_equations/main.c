#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <time.h>

double randomNumber( unsigned int *seed ){     //generate random numbers to test program.
	double maxRand = (double)RAND_MAX;
	double randNum = (double)rand_r( seed );
	return randNum/maxRand;
}

void set_data_matrix(gsl_matrix* dat_matrix, gsl_vector* dat_vector, unsigned int *seed){ //make matrix
	for(int row_i = 0; row_i < dat_matrix -> nrow; row_i++){
		for(int col_i = 0; col_i < dat_matrix -> ncol; col_i++){
			gsl_matrix_set(dat_matrix, row_i, col_i, randomNumber(seed));
		}
		gsl_vector_set(dat_vector, row_i, randomNumber(seed));
	}
}




void GS_decomp(gsl_matrix* A, gsl_matrix* R){
// GS orth. on a matrix (A), returns Q and R.

	gsl_vector* a_i = gsl_vector_alloc(A->nrow);
	gsl_vector* q_i = gsl_vector_alloc(A->nrow);


	for (int i=0; i<(A->ncol); i++){ // run for each col
		gsl_matrix_get_col(a_i, A, i);
		gsl_vector_memcpy(q_i,a_i);


		double a_i_norm = gsl_blas_dnrm2(a_i); //compute norm
		gsl_matrix_set(R,i,i, a_i_norm);
		gsl_vector_scale(q_i,1/a_i_norm);

		gsl_matrix_set_col(A,i,q_i);

		for (int j=i+1; j<(A->ncol); j++){ // loop over col right of i
			gsl_matrix_get_col(a_i,A,j);
			double R_ij=0; //define [i,j] element in matrix
			gsl_blas_ddot(q_i,a_i,&R_ij);
			gsl_matrix_set(R,i,j,R_ij);
			gsl_matrix_set(R,j,i,0.);
			gsl_blas_daxpy(-R_ij,q_i,a_i);
			gsl_matrix_set_col(A,j,a_i);
		}
	}

		gsl_vector_free(a_i);
		gsl_vector_free(q_i);
}




int main(){




return 0;
}






