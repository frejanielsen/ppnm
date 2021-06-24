#ifndef HAVE_RUNGEKUTTA_H
#define HAVE_RUNGEKUTTA_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void
rk_stepXY(void (*f)(double t, gsl_vector* y, gsl_vector* dydt), double var, 
	gsl_vector *f_val, double step, gsl_vector* f_step, gsl_vector *err);

void rk_driver(void (*f)(double t, gsl_vector * y, gsl_vector * dydt), double a,
	gsl_vector * f_val_a, double b, gsl_vector *f_val_b, double step,
	double abs_acc, double rel_acc, FILE *path_to_file);

#endif //HAVE_RUNGEKUTTA_H
