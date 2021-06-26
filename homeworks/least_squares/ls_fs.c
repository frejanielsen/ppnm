#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>


void matrixprint(gsl_matrix* A){
	int n=A->size1;
	int m=A->size2;

	for(int i=0;i<n;i++){
	        for(int j=0;j<m;j++){
	                printf("%10.2g",gsl_matrix_get(A,i,j));
	        }
	        printf("\n");
}

printf("\n\n");
}


void* GS_decomp(gsl_matrix* A, gsl_matrix* R){
	//replace A bu Q,R
	int size_c=A->size2; //number of columns in A
	for (int i=0;i<size_c;i++){
		gsl_vector_view ai = gsl_matrix_column(A,i);
		double ai_norm= gsl_blas_dnrm2(&ai.vector);
		gsl_matrix_set(R,i,i,ai_norm);
		gsl_vector_scale(&ai.vector,1/ai_norm);

		for (int j=i+1;j<size_c;j++){
			gsl_vector_view aj = gsl_matrix_column(A,j);
			double proj; //holder for projection value
			gsl_blas_ddot(&ai.vector,&aj.vector,&proj);
			gsl_blas_daxpy(-proj,&ai.vector,&aj.vector);
			gsl_matrix_set(R,i,j,proj); //inserting valuje in upper triangle of R
			gsl_matrix_set(R,j,i,0); //setting entrance in lower triangle to 0
		}

	}


return NULL;

}


void*  backsub(gsl_matrix* R, gsl_vector* x){
	int size_row=R->size1;
	for (int i=size_row-1; i>=0; i--){
		double ph=0; //holder for value for the backsubstitution
		for(int j=i+1; j<size_row;j++){
			ph+=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		}
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-ph)/gsl_matrix_get(R,i,i));
	}

return NULL;
}

void* GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
	backsub(R,x);
return NULL;
}


void* invert(gsl_matrix* Q,gsl_matrix* R,gsl_matrix* Ai){
        int n=Q->size1;
        gsl_vector* ei=gsl_vector_alloc(n);
        for (int i=0;i<n;i++){
                gsl_vector_set(ei,i,0);
        }
        gsl_vector_set(ei,0,1);
        gsl_vector* aii=gsl_vector_alloc(n);
        for (int i=0;i<n;i++){
                gsl_vector_set(aii,i,gsl_matrix_get(Ai,i,0));
        }
        GS_solve(Q,R,ei,aii);

        for (int i=0;i<n;i++){
                gsl_matrix_set(Ai,i,0,gsl_vector_get(aii,i));
        }
        for(int i=1;i<n;i++){
                gsl_vector_set(ei,i,1);
                gsl_vector_set(ei,i-1,0);
                for (int i=0;i<n;i++){
                        gsl_vector_set(aii,i,gsl_matrix_get(Ai,i,0));
                }
                GS_solve(Q,R,ei,aii);
                for (int j=0;j<n;j++){
                        gsl_matrix_set(Ai,j,i,gsl_vector_get(aii,j));
                }
	}
	gsl_vector_free(aii);
	return NULL;
}



double func(int i, double x){
	switch(i){
		case  0: return 1  ; break;
		case  1: return -x  ; break;
		case  2: return x*x; break;
		default: return NAN;
		}
	}

void* fit(gsl_matrix* data, double func(int i, double x), int N_f, gsl_vector* y, gsl_matrix* cov ){
        int N = data->size1;
	gsl_matrix* A = gsl_matrix_alloc(N,N_f);
	gsl_matrix* R = gsl_matrix_alloc(N_f,N_f);
	gsl_matrix* Rr = gsl_matrix_alloc(N_f,N_f);
	gsl_matrix* Ri = gsl_matrix_alloc(N_f,N_f);
	gsl_vector* Y = gsl_vector_alloc(N);

	for (int i=0;i<N;i++){
		for (int j=0;j<N_f;j++){
			gsl_matrix_set(A,i,j,func(j,gsl_matrix_get(data,i,0))/gsl_matrix_get(data,i,2));
		}
	}
	matrixprint(A);
	for (int j=0;j<N;j++){
        	gsl_vector_set(Y,j,gsl_matrix_get(data,j,1)/gsl_matrix_get(data,j,2));
        }

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,R);
	GS_decomp(R,Rr);
	invert(R,Rr,cov);

	GS_decomp(A,R);
        gsl_blas_dgemv(CblasTrans,1.0,A,Y,0.0,y);
        backsub(R,y);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(Rr);
	gsl_matrix_free(Ri);
	gsl_vector_free(Y);


return NULL;
}


