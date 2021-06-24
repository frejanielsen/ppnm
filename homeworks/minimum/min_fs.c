#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 3e-16
#endif
static const double DELTA=sqrt(DBL_EPSILON);

void gradient(double F(gsl_vector* x),gsl_vector* x, gsl_vector* dF){

int n=x->size;
double Fx=F(x);
double Fx_dx;
double dx;
double xi;
for (int i=0;i<n;i++){
	xi=gsl_vector_get(x,i);
	if (fabs(xi)<sqrt(DELTA)){
		dx=DELTA;
	} else {
		dx=fabs(xi)*DELTA;
	}
        gsl_vector_set(x,i,gsl_vector_get(x,i)+dx);
        Fx_dx = F(x);
        gsl_vector_set(x,i,gsl_vector_get(x,i)-dx);
        gsl_vector_set(dF,i,(Fx_dx-Fx)/dx);

	}
}

void qnewton(double F(gsl_vector* x), gsl_vector* x, double eps){

int dim = x->size;
int nsteps = 0;
int ngood = 0;
int nbad = 0;

gsl_matrix* B = gsl_matrix_alloc(dim,dim);
gsl_vector* dF = gsl_vector_alloc(dim);
gsl_vector* dx = gsl_vector_alloc(dim);
gsl_vector* x_dx = gsl_vector_alloc(dim);
gsl_vector* dx_dx = gsl_vector_alloc(dim);
gsl_vector* y = gsl_vector_alloc(dim);
gsl_vector* u = gsl_vector_alloc(dim);
//gsl_vector* a = gsl_vector_alloc(dim);
gsl_matrix_set_identity(B);

gradient(F,x,dF);
double Fx = F(x);
double Fx_dx;

while (nsteps<100){
	nsteps++;
	gsl_blas_dgemv(CblasNoTrans,-1.0,B,dF,0.0,dx);
	if (gsl_blas_dnrm2(dx)<DELTA*gsl_blas_dnrm2(x) || gsl_blas_dnrm2(dF)<eps){
		break;
	}
	double lambda = 1;
	while(1){
		gsl_vector_memcpy(x_dx,x);
		gsl_vector_add(x_dx,dx);
		Fx_dx = F(x_dx);
		double sdotdF;
		gsl_blas_ddot(dx,dF,&sdotdF);
		if (Fx_dx<Fx+0.01*sdotdF){
			ngood++; break;}
		if (lambda<DELTA){
			nbad++; gsl_matrix_set_identity(B); break; }
		lambda /= 2;
		gsl_vector_scale(dx,0.5);
	}

	gradient(F, x_dx, dx_dx);
	gsl_vector_memcpy(y,dx_dx);
	gsl_blas_daxpy(-1,dF,y); //sum
	gsl_vector_memcpy(u,dx);
	gsl_blas_dgemv(CblasNoTrans,1.0,B,y,1.0,u);
	double sdoty;
	double udoty;
	gsl_blas_ddot(dx,y,&udoty);
	if (fabs(sdoty)>1e-12){
		gsl_blas_ddot(u,y,&udoty);
		double gamma = udoty/2/sdoty;
		gsl_blas_daxpy(-1.0*gamma,dx,u);
		gsl_blas_dger(1/sdoty,u,dx,B);
		gsl_blas_dger(1/sdoty,dx,u,B);
	}
	gsl_vector_memcpy(x,x_dx);
	gsl_vector_memcpy(dF,dx_dx);
	Fx = Fx_dx;
}

gsl_matrix_free(B);
gsl_vector_free(dF);
gsl_vector_free(dx);
gsl_vector_free(x_dx);
gsl_vector_free(dx_dx);
gsl_vector_free(y);
gsl_vector_free(u);
//gsl_vector_free(a);

}

















