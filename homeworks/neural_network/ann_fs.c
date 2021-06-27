
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX;
/*
typedef struct { int n; double(*f)(double); gsl_vector* params; } ann;

ann* alloc(int n,double(*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n = n;
	network->f = f;
	network->params = gsl_vector_alloc(3*n);
	return network;
}


void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}


double ann_response(ann* network, double x){
	int N = network->n;
	double sum = 0;
	for(int i=0;i<N;i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i + 1);
		double w = gsl_vector_get(network->params,3*i + 2);
		sum += network->f((x-a)/b)*w;
	}
	return sum;
}

ann* network;


double cost(gsl_vector* p){
	ann* network; gsl_vector* xs; gsl_vector*ys;
	gsl_vector_memcpy(network->params,p);
	int N = xs->size;;
	double sum = 0;
	for(int i =0;i<N;i++){
		double xi = gsl_vector_get(xs,i);
		double yi = gsl_vector_get(ys,i);
		double fi = ann_response(network,xi);
		sum += fabs(fi-yi)*fabs(fi- yi);
	}
	return sum / N;
}

*/
void gradient(double f(gsl_vector* xvector),gsl_vector* x,gsl_vector* grad){
	double DELTA = 2.22045e-16;
	int n = x->size;
	double dx = sqrt(DELTA);
	double fdx = 0;
	double fx = 0;
	double xi = 0;
	double gradi =0;
	for(int i=0;i<n;i++){
		xi = gsl_vector_get(x,i);
		fx = f(x);
		gsl_vector_set(x,i,xi+dx);
		fdx = f(x);
		gradi = (fdx - fx ) / ( dx );
		gsl_vector_set(grad,i,gradi);
		gsl_vector_set(x,i,xi);
		}
}



void multimc(double f(gsl_vector* xvector),gsl_vector* x, double eps,int steps){
	int n = x->size;
	double DELTA = 2.22045e-16;
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* grad = gsl_vector_alloc(n);
	gsl_vector* gradxs = gsl_vector_alloc(n);
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_vector* xs = gsl_vector_alloc(n);
	gsl_vector* s = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_matrix_set_identity(B);
	gsl_matrix_set_identity(I);
	double fx = f(x);
	double fxs = 0;
	double sTg=0;
	gradient(f,x,grad);
	int k = 0;
	while(k<steps){
		k++;
		gsl_blas_dgemv(CblasNoTrans, -1, B, grad,0,Dx); // eq(6)
		double lambda = 1;
		double norm_grad = gsl_blas_dnrm2(grad);
		if(norm_grad<DELTA*gsl_blas_dnrm2(x)){
			//printf("The gradient is within our acceptance\nit is %g",norm_grad);
			break;
		}
		double norm_Dx = gsl_blas_dnrm2(Dx);
		if(norm_Dx<eps){
			//printf("Dx is within our acceptance\n");
			break;
		}

		while(1){

			gsl_vector_memcpy(s,Dx);
			gsl_vector_scale(s,lambda);
			if(gsl_blas_dnrm2(s)==0){
				printf("s is zero\n");
				break;
			}
			gsl_vector_memcpy(xs,x);
			gsl_vector_add(xs,s); // x+ s
			fxs = f(xs); // f(x+s)
			gsl_blas_ddot(s,grad,&sTg);
			if(fxs<fx+0.001*sTg){
				break;
				}
			if(lambda<0.01){
				gsl_matrix_set_identity(B);
				break;
			}
			lambda *=0.5;
			}
			if(gsl_blas_dnrm2(Dx)<1e-10){
				break;
			}
		gradient(f,xs,gradxs);
		gsl_vector_memcpy(y,gradxs);
		gsl_vector_memcpy(u,s);
		gsl_vector_sub(y,grad); // y = Grad(x+s) - Grad (x)
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u);
		double sTy = 0;
		gsl_blas_ddot(s,y,&sTy);
		if(fabs(sTy)>10e-6){
			double gamma = 0;
			double uTy = 0;
			gsl_blas_ddot(u,y,&uTy); // uTy
			gamma = 0.5 *(uTy)*(1/sTy); // sTy
			gsl_vector_memcpy(c,u);
			gsl_blas_daxpy(-gamma,s,u); // sTy * a
			gsl_vector_scale(u,1/sTy);
			gsl_blas_dger(1,u,s,B);
			gsl_blas_dger(1,s,u,B);
		}
		gsl_vector_memcpy(x,xs);
		gsl_vector_memcpy(grad,gradxs);
		fx = fxs;
	}
	//printf("The number of steps it took to find the extremum is %i\n",k);
	gsl_matrix_free(I);
	gsl_matrix_free(B);
	gsl_vector_free(Dx);
	gsl_vector_free(grad);
	gsl_vector_free(gradxs);
	gsl_vector_free(xs);
	gsl_vector_free(s);
	gsl_vector_free(c);
	gsl_vector_free(u);
	gsl_vector_free(y);
}

/*
void ann_train(ann* network,gsl_vector* xs, gsl_vector* ys){
	int N = network->params->size;
	gsl_vector* p = gsl_vector_alloc(N);
	gsl_vector_memcpy(p,network->params);

	multimc(cost,p,1e-3,1000);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}

*/
