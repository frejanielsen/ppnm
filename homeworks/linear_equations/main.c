#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include <assert.h>

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
	int n=A->size1;
		int m=A->size2;

		for(int i=0;i<n;i++){
		        for(int j=0;j<m;j++){
		                fprintf(file,"%10.2g",gsl_matrix_get(A,i,j));
		        }
		        fprintf(file,"\n");}
}


void print_vector(gsl_vector* v, FILE* file){
	for(int i=0;i<v->size;i++)fprintf(file, "%10g ",gsl_vector_get(v,i));
}




void GS_decomp(gsl_matrix* A, gsl_matrix* R){
// GS orth. on a matrix (A), returns Q and R.
/*
	gsl_vector* a_i = gsl_vector_alloc(A->size1);
	gsl_vector* q_i = gsl_vector_alloc(A->size2);
*/
	int size_c=A->size2; //number of columns in A

	for (int i=0;i<size_c;i++){
		gsl_vector_view a_i = gsl_matrix_column(A,i);
		double a_i_norm= gsl_blas_dnrm2(&a_i.vector);
		gsl_matrix_set(R,i,i,a_i_norm);
		gsl_vector_scale(&a_i.vector,1/a_i_norm);

			for (int j=i+1;j<size_c;j++){
				gsl_vector_view a_j = gsl_matrix_column(A,j);
				double proj; //holder for projection value
				gsl_blas_ddot(&a_i.vector,&a_j.vector,&proj);
				gsl_blas_daxpy(-proj,&a_i.vector,&a_j.vector);
				gsl_matrix_set(R,i,j,proj); //inserting valuje in upper triangle of R
				gsl_matrix_set(R,j,i,0); //setting entrance in lower triangle to 0
			}

	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){

	gsl_blas_dgemv(CblasTrans,1.,Q,b,0.,x);  // matrix product of Q^T and b  saves in x

	//perform back substirution on R - return in x - perhaps in own function
	int size_row=R->size1;
		for (int i=size_row-1; i>=0; i--){
			double ph=0; //holder for value for the backsubstitution
			for(int j=i+1; j<size_row;j++){
				ph+=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
			}
			gsl_vector_set(x,i,(gsl_vector_get(x,i)-ph)/gsl_matrix_get(R,i,i));
		}

}


void* GS_inverse(gsl_matrix* Q, gsl_matrix* R,gsl_matrix* B){
	int n=Q->size1;
	gsl_vector* e = gsl_vector_alloc(n);
	for (int i=0;i<n;i++){
		gsl_vector_set(e,i,0);
	}
	gsl_vector_set(e,0,1);
	gsl_vector* bii = gsl_vector_alloc(n);
	for (int i=0;i<n;i++){
		gsl_vector_set(bii, i, gsl_matrix_get(B,i,0));
	}
	GS_solve(Q, R, e, bii);

	for (int i=0;i<n;i++){
		gsl_matrix_set(B,i,0,gsl_vector_get(bii,i));
	}
	for(int i=1;i<n;i++){
		gsl_vector_set(e,i,1);
		gsl_vector_set(e,i-1,0);
		for (int i=0;i<n;i++){
			gsl_vector_set(bii,i,gsl_matrix_get(B,i,0));
		}
		GS_solve(Q, R, e, bii);
		for (int j=0;j<n;j++){
			gsl_matrix_set(B,j,i,gsl_vector_get(bii,j));
		}
	}

gsl_vector_free(bii);
return NULL;
}


int main(){
	// A.1
	int size[2];

	size[0]=7;
	size[1]=5;
	assert(size[0]>=size[1]);

	FILE* outA = fopen("out_lin.txt","w");

	// allocate menory for matricies
	gsl_matrix* A = gsl_matrix_alloc(size[0],size[1]);
	gsl_matrix* R = gsl_matrix_alloc(size[1],size[1]);
	gsl_matrix* A_c = gsl_matrix_alloc(size[0],size[1]);

	/*for(int i=0; i<(A->size1);i++){ // random [n,m] matrix
		for(int j=0; j<(A->size2); j++){
			double Aij =rand()/100;
			gsl_matrix_set(A,i,j,Aij);
		}
	}*/

	unsigned int rand_seed=2;
	for (int i = 0; i<size[1]; i++){
		for (int j=0; j<size[0]; j++){
			gsl_matrix_set(A,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		}}

	//gsl_matrix_memcpy(A_c,A); // copy A

	fprintf(outA,"Part A.1:\n \n");
	fprintf(outA," A=\n");
	print_matrix(A,outA);

	fprintf(outA,"QR decomposition\n");
	GS_decomp(A,R); //GS on A

	fprintf(outA,"\n R is upper triangular? \n");
	print_matrix(R,outA);

	gsl_matrix* Qt = gsl_matrix_alloc(size[1],size[1]); // allocate menory for QT*Q
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,A,A,0,Qt);  //Q^T*Q
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A,R,-1.,A_c);

	fprintf(outA," Q=\n");
	print_matrix(A,outA);


	fprintf(outA,"\n is Q^T*Q=1 ?\n");
	print_matrix(Qt,outA);
	fprintf(outA,"\n is QR=A i.e. QR-A=0 ?\n");
	print_matrix(A_c,outA);
	fprintf(outA,"\n \n \n");

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(A_c);
	//A.2
	A = gsl_matrix_alloc(size[1],size[1]);
	R = gsl_matrix_alloc(size[1],size[1]);
	A_c = gsl_matrix_alloc(size[1],size[1]);
	gsl_vector* b = gsl_vector_alloc(size[1]);
	gsl_vector* x = gsl_vector_alloc(size[1]);

	for (int i = 0; i<size[1]; i++){
		for (int j=0; j<size[1]; j++){
			gsl_matrix_set(A,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		}
		gsl_vector_set(b,i,1.0*rand_r(&rand_seed)/RAND_MAX);
	}


	/*
	for(int i=0; i<(b->size); i++){
		double bi = rand()/100;
		gsl_vector_set(b,i,bi);
	}
	*/
	fprintf(outA,"Part A.2:\n \n");
	fprintf(outA," A=\n");
	print_matrix(A,outA);
	fprintf(outA," b=\n");
	print_vector(b,outA);

	//fprintf(outA,"QR decomposition\n");

	gsl_matrix_memcpy(A_c,A); // copy A
	GS_decomp(A,R); //GS on A
	GS_solve(A,R,b,x); //slves O*R*x=b


	fprintf(outA, "\n Is A*x = b? A*x=\n");
	gsl_blas_dgemv(CblasNoTrans,1.,A_c,x,0.0,b);
	print_vector(b,outA);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(A_c);

	//B
	//FILE* outB = fopen("out_lin.txt","w");
	fprintf(outA,"\n \n");
	fprintf(outA,"Part B:\n \n");
	int n = 6;

	A = gsl_matrix_alloc(n,n);
	R = gsl_matrix_alloc(n,n);
	A_c = gsl_matrix_alloc(n,n);

	for (int i = 0; i<n; i++){
		for (int j=0; j<n; j++){
			gsl_matrix_set(A,j,i,1.0*rand_r(&rand_seed)/RAND_MAX);
		}
	}


	fprintf(outA," A=\n");
	print_matrix(A,outA);

	gsl_matrix_memcpy(A_c,A); // copy A
	GS_decomp(A,R);

	gsl_matrix* A_in = gsl_matrix_alloc(n,n);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	GS_inverse(A,R,A_in);


	fprintf(outA," A^-1=B=\n");
	print_matrix(A_in,outA);


	fprintf(outA,"check that BA=I=AB\n B*A=\n ");

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_in,A_c,0.0,I);

	for (int i =0; i<n;i++){
		for(int j=0; j<n;j++){
			fprintf(outA,"%10g   ",gsl_matrix_get(I,i,j));
		}
		fprintf(outA,"\n");
	}


	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_c,A_in,0.0,I);
	fprintf(outA,"A*B=\n");
	for (int i =0; i<n;i++){
		for(int j=0; j<n;j++){
			fprintf(outA,"%10g   ",gsl_matrix_get(I,i,j));
		}
		fprintf(outA,"\n");
	}


	//C:


	clock_t t = clock();

	FILE* my_time = fopen("out_my_time.txt","w");
	FILE* gsl_time = fopen("out_gsl_time.txt","w");
	double time_taken;

	int nmin=300;
	int nmax=800;
	int nstep=10;

	for(int n=nmin;n<=nmax;n+=nstep){
		gsl_matrix* A = gsl_matrix_alloc(n,n);
		gsl_matrix* R = gsl_matrix_alloc(n,n);
		for (int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				gsl_matrix_set(A,i,j,rand_r(&rand_seed)/RAND_MAX);
			 }
		}

		t = clock();
		GS_decomp(A,R);
		t = clock()-t;
		time_taken = ((double) t)/CLOCKS_PER_SEC;
		fprintf(my_time,"%10i  %10g\n",n,time_taken);
		gsl_matrix_free(A);
		gsl_matrix_free(R);
	}

	nmin=2000;
	nmax=4000;
	nstep=50;

	for(int n=nmin;n<=nmax;n+=nstep){
		gsl_matrix* A = gsl_matrix_alloc(n,n);
		gsl_vector* R = gsl_vector_alloc(n);
		for (int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				gsl_matrix_set(A,i,j,rand_r(&rand_seed)/RAND_MAX);
			 }
		}

		t=clock();
		gsl_linalg_QR_decomp(A,R);
		t=clock()-t;
		time_taken=((double) t)/CLOCKS_PER_SEC;
		fprintf(gsl_time,"%10i  %10g\n",n,time_taken);
		gsl_matrix_free(A);
		gsl_vector_free(R);
	}






return 0;
}

