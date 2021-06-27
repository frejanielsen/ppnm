#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


#include "ann_fs.h"
#include "printfs.h"



double func_act(double x){
	return x * exp(-x*x);
}

double func(double x){
	return sin(5*x+2)*exp(-x*x);
}
/*
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
*/

ann* network;


int main(){
	int N = 6;
	int M = 3*N;
	//ann* network;
	gsl_vector* xs = gsl_vector_alloc(M);
	gsl_vector* ys = gsl_vector_alloc(M);
	network = ann_alloc(N,func_act);
	double a0 = -1;
	double an = 1;


	for(int i = 0; i<M;i++){
		double xi = a0 + (an-a0)*i/(M-1);
		gsl_vector_set(xs,i,xi);
		double fi = func(xi);
		gsl_vector_set(ys,i,fi);
		fprintf(stderr,"%g %g \n",xi,fi);
	}

	gsl_vector_set(network->params,3*N-1,1);
	for(int i = 0;i<network->n;i++){
		gsl_vector_set(network->params,1 * i ,a0 + (an -a0)*i/(network->n-1));
		gsl_vector_set(network->params,3 * i + 1,i);
		gsl_vector_set(network->params,3 * i + 2,i);
	}

	ann_train(network,xs,ys);
	FILE* data = fopen("output.txt","w");
	fprintf(data,"data:\n");
	for(double i =a0;i<an;i+=1.0/72){
		double expct_output = ann_response(network, i);
		fprintf(data,"%g %g \n",i,expct_output);
	}

	fprintf(data,"\n\n original\n");
	for(double i =0;i<M;i++){
		fprintf(data,"%g %g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i));
	}
	fclose(data);
	ann_free(network);
	gsl_vector_free(xs);
	gsl_vector_free(ys);













return 0;
}
