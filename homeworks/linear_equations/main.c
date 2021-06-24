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
	for(int row_i = 0; row_i < dat_matrix -> size1; row_i++){
		for(int col_i = 0; col_i < dat_matrix -> size2; col_i++){
			gsl_matrix_set(dat_matrix, row_i, col_i, randomNumber(seed));
		}
		gsl_vector_set(dat_vector, row_i, randomNumber(seed));
	}
}

void print_matrix(gsl_matrix* A,FILE* file){
	for(int i=0; i<(A->size1);i++){
		for(int j=0; j<(A->size2); j++){
			double Aij = gsl_matrix_get(A,i,j);
			fprintf(file, "%0.3g ",Aij);
		}
		fprintf(file,"\n");
	}
}

void print_vector(gsl_vector* v, FILE* file){
	for(int i=0;i<v->size;i++)fprintf(file, "%10g ",gsl_vector_get(v,i));
}




void GS_decomp(gsl_matrix* A, gsl_matrix* R){
// GS orth. on a matrix (A), returns Q and R.

	gsl_vector* a_i = gsl_vector_alloc(A->size1);
	gsl_vector* q_i = gsl_vector_alloc(A->size2);


	for (int i=0; i<(A->size2); i++){ // run for each col
		gsl_matrix_get_col(a_i, A, i);
		gsl_vector_memcpy(q_i,a_i);


		double a_i_norm = gsl_blas_dnrm2(a_i); //compute norm
		gsl_matrix_set(R,i,i, a_i_norm);
		gsl_vector_scale(q_i,1/a_i_norm);

		gsl_matrix_set_col(A,i,q_i);

		for (int j=i+1; j<(A->size2); j++){ // loop over col right of i
			gsl_matrix_get_col(a_i,A,j);
			double R_ij=0; //define [i,j] element in matrix
			gsl_blas_ddot(q_i,a_i,&R_ij);
			gsl_matrix_set(R,i,j,R_ij);
			gsl_matrix_set(R,j,i,0.);    // make R upper triangular
			gsl_blas_daxpy(-R_ij,q_i,a_i);
			gsl_matrix_set_col(A,j,a_i);
		}
	}

		gsl_vector_free(a_i);
		gsl_vector_free(q_i);
}



void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){

	gsl_blas_dgemv(CblasTrans,1.,Q,b,0.,x);  // matrix product of Q^T and b  saves in x

	//perform back substirution on R - return in x - perhaps in own function
	for(int i=(x->size)-1; i>=0; i--){
		double R_ii = gsl_matrix_get(R,i,i);
		double x_i = gsl_vector_get(x,i);
		for (int j=i+1; j<(x->size);j++){
			double R_ij = gsl_matrix_get(R,i,j);
			double x_j = gsl_vector_get(x,j);
			x_i -= R_ij*x_j;
		}

		x_i /= R_ii;
		gsl_vector_set(x,i,x_i);
	}


}


int main(){

	// A.1

	int n = 8;
	int m = 8;

	FILE* outA = fopen("out.A.txt","w");

	// allocate menory for matricies
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* A_c = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);

	for(int i=0; i<(A->size1);i++){ // random [n,m] matrix
		for(int j=0; j<(A->size2); j++){
			double Aij =rand()/100;
			gsl_matrix_set(A,i,j,Aij);
		}
	}

	gsl_matrix_memcpy(A_c,A); // copy A
	GS_decomp(A,R); //GS on A

	gsl_matrix* Qt = gsl_matrix_alloc(m,m); // allocate menory for QT*Q
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,A,A,0,Qt);  //Q^T*Q
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A,R,-1.,A_c);

	fprintf(outA," R is upper triangular? \n");
	print_matrix(R,outA);
	fprintf(outA,"\n is Q^T*Q=1 ?\n");
	print_matrix(Qt,outA);
	fprintf(outA,"\n is QR=A i.e. QR-A=0 ?\n");
	print_matrix(A_c,outA);


	//A.2
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(m);


	for(int i=0; i<(A->size1);i++){ // random [n,m] matrix
		for(int j=0; j<(A->size2); j++){
			double Aij =rand()/100;
			gsl_matrix_set(A,i,j,Aij);
		}
	}

	for(int i=0; i<(b->size); i++){
		double bi = rand()/100;
		gsl_vector_set(b,i,bi);
	}

	gsl_matrix_memcpy(A_c,A); // copy A
	GS_decomp(A,R); //GS on A
	GS_solve(A,R,b,x); //slves O*R*x=b
	gsl_blas_dgemv(CblasNoTrans,1.,A_c,x,-1.,b); //A*x-b

	fprintf(outA, "\n Is A*x = b? ie. A*x-b=0 \n");
	print_vector(b,outA);


	//B
	FILE* outB = fopen("out.B.txt","w");


return 0;
}






