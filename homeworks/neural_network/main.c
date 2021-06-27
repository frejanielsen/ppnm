#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>

#include "ann_fs.h"



double func_act(double x){
	return x * exp(-x*x);
}

double func(double x){
	return sin(5*x+2)*exp(-x*x);
}

typedef struct{
	int n;
       	double (*f)(double);
       	gsl_vector* params;
	} ann;

ann* ann_alloc(int n, double (*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network ->f=f;
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

ann* network; gsl_vector* xs; gsl_vector*ys;



double cost(gsl_vector* p){
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


void ann_train(ann* network,gsl_vector* xs, gsl_vector* ys){
	int N = network->params->size;
	gsl_vector* p = gsl_vector_alloc(N);
	gsl_vector_memcpy(p,network->params);

	multimc(cost,p,1e-3,1000);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}


int main(){
	int N = 6;
	int M = 3*N;
	//ann* network;
	gsl_vector* xs = gsl_vector_alloc(M);
	gsl_vector* ys = gsl_vector_alloc(M);
	network = ann_alloc(N,func_act);
	double a0 = -1;
	double an = 1;

	// generate x and y
	for(int i = 0; i<M;i++){
		double xi = a0 + (an-a0)*i/(M-1);
		gsl_vector_set(xs,i,xi);
		double fi = func(xi);
		gsl_vector_set(ys,i,fi);
		//fprintf(stderr,"%g %g \n",xi,fi);
	}
	//initial params
	gsl_vector_set(network->params,3*N-1,1);
	for(int j = 0;j<network->n;j++){
		gsl_vector_set(network->params,1 * j ,a0 + (an -a0)*j/(network->n-1));
		gsl_vector_set(network->params,3 * j + 1,j);
		gsl_vector_set(network->params,3 * j + 2,j);
	}

	ann_train(network,xs,ys);
	FILE* out = fopen("output.txt","w");

	for(int i=0; i<M; i++) {
			double x = gsl_vector_get(xs, i);
			double f = gsl_vector_get(ys, i);
			double val = ann_response(network,x);
			fprintf(out,"x=%10g, f(x)=%10g, NN(x)=%10g\n", x, f, val);
		}
		fprintf(out,"\n\n");



/*
	FILE* data = fopen("output_ann.txt","w");

	fprintf(data,"#data:\n");

	for(double i = a0; i < an; i+=1.0/72){
		double expct_output = ann_response(network, i);
		fprintf(data,"%g %g \n",i,expct_output);
	}

	fprintf(data,"\n\n # original:\n");
	for(double k =0;k<M;k++){
		fprintf(data,"%g %g\n",gsl_vector_get(xs,k),gsl_vector_get(ys,k));
	}
	fclose(data);
*/
	ann_free(network);
	gsl_vector_free(xs);
	gsl_vector_free(ys);
return 0;
}
