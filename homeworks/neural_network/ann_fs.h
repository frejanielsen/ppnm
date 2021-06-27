#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_linalg.h>
#include<stdio.h>
#define rnd (double)rand()/RAND_MAX
#ifndef HAVE_ANN_H
#define HAVE_ANN_H
typedef struct { int n; double(*f)(double); gsl_vector* params; } ann;
ann* ann_alloc(int n, double (*f)(double));
void ann_free(ann* network);
double ann_response(ann* network, double x);
double cost(gsl_vector* p);
void gradient(double f(gsl_vector* xvector),gsl_vector* x,gsl_vector* grad);
void multimc(double f(gsl_vector* xvector),gsl_vector* x, double eps,int steps);
void ann_train(ann* list, gsl_vector* xs, gsl_vector* ys);

#endif
