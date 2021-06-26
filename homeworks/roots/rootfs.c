#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#include "linfs.h"

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){

double d_x = sqrt(3e-16); //0.0001;

int dim = x->size;
double lambda = 1;
double lambda_min = 0.00001;
int c = 0;


gsl_vector* fdx = gsl_vector_alloc(dim);
gsl_vector* fx = gsl_vector_alloc(dim);
gsl_vector* dx = gsl_vector_alloc(dim);
gsl_matrix* J = gsl_matrix_alloc(dim,dim);
gsl_matrix* R = gsl_matrix_alloc(dim,dim);

for (int i = 0; i < dim; i++){
	gsl_vector_set(dx,i,0);
}

do{
	c++;
	//assert(c<10000);
	lambda=1;
	f(x,fx);
	gsl_blas_dscal(0,dx);
	gsl_vector_set(dx,0,d_x);
	gsl_blas_daxpy(1.0,dx,x);
	f(x,fdx);
	gsl_blas_daxpy(-1.0,dx,x);
	for (int i=0; i<dim;i++){
		gsl_matrix_set(J,i,0,(gsl_vector_get(fdx,i)-gsl_vector_get(fx,i))/d_x);
	}

	for (int i=1;i<dim;i++){
		gsl_vector_set(dx,i-1,0);
		gsl_vector_set(dx,i,d_x);

		gsl_blas_daxpy(1.0,dx,x);
		f(x,fdx);
		gsl_blas_daxpy(-1.0,dx,x);
		for (int j=0; j<dim;j++){
        	        gsl_matrix_set(J,j,i,(gsl_vector_get(fdx,j)-gsl_vector_get(fx,j))*1.0/d_x);
	        }
	}

	GS_decomp(J,R);
	gsl_blas_daxpy(-2.0,fx,fx);
	GS_solve(J,R,fx,dx);
	gsl_blas_daxpy(-2.0,fx,fx);
	gsl_blas_daxpy(1.0,dx,x);
	f(x,fdx);
	gsl_blas_daxpy(-1.0,dx,x);
	while(gsl_blas_dnrm2(fdx)>(1-lambda/2)*gsl_blas_dnrm2(fx) && lambda>lambda_min){
		lambda/=2;
	}
	gsl_blas_daxpy(lambda,dx,x);
	f(x,fx);

} while (gsl_blas_dnrm2(fx)>eps);

gsl_vector_free(dx);
gsl_vector_free(fx);
gsl_vector_free(fdx);
gsl_matrix_free(J);
gsl_matrix_free(R);

}

