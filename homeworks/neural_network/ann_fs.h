#ifndef HAVE_NEURALNETWORK_H
#define HAVE_NEURALNETWORK_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double random_number(unsigned int *seed);

void vectorprint(char *s, gsl_vector *v);

void matrixprint(int n, gsl_matrix *A, char *s);

void symmetric(gsl_matrix *A, unsigned int *seed);


void gradient(double f(gsl_vector *), gsl_vector *m, gsl_vector *g);

void multimc(double f(gsl_vector *), gsl_vector *m, double tol);



typedef struct
{
    int n;

    double (*f)(double);

    double (*df)(double);

    double (*fi)(double);

    gsl_vector *params;
} ann;

double ann_response(ann *network, double ep);

void ann_train(ann *network, gsl_vector *data, gsl_vector *labels);

double ann_response_d(ann *network, double ep);

double ann_response_i(ann *network, double rp, double lp);

ann*ann_alloc(int n, double (*f)(double), double (*df)(double), double (*fi)(double));

void ann_free(ann *network);


#endif //HAVE_NEURALNETWORK_H

