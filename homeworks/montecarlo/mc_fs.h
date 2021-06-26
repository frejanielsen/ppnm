#include <complex.h>
//complex double plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N);
void plainmc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* err);
void multimc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* eps);
double corput(int,int);
