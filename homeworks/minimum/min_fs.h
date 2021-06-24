#ifndef MIN_FS
#define MIN_FS
void gradient(double F(gsl_vector* x),gsl_vector* x, gsl_vector* dF);
void qnewton(double F(gsl_vector* x), gsl_vector* x, double eps);

#endif
