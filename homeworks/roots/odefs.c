#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
//#include "rung_kutta.h"
//#include "diff_eqs.h"



void rk_stepXY(void f (double t,gsl_vector* y,gsl_vector* dydt),
	double var, gsl_vector* f_val, double step, gsl_vector* f_step,
	gsl_vector* err){

	int order = f_val->size; //Order of differential equation
	gsl_vector *tan0 = gsl_vector_alloc(order); //k_0
	gsl_vector *tan12 = gsl_vector_alloc(order); //k_1/2
	gsl_vector *temp_f_val = gsl_vector_alloc(order);

	f(var, f_val, tan0); //RHS of ODE

	for (int i = 0; i < order; i++){
		double f_val_i = gsl_vector_get(f_val, i);
		double tan0i = gsl_vector_get(tan0, i);
		double temp_f_val_i = f_val_i + 0.5 * tan0i * step;
		gsl_vector_set(temp_f_val, i, temp_f_val_i);
	}

	f(var + step * 0.5, temp_f_val, tan12);

	//Advance the solution
	for (int i = 0; i < order; i++){
		double f_val_i = gsl_vector_get(f_val, i);
		double tan12i = gsl_vector_get(tan12, i);
		double temp_f_step_i = f_val_i + tan12i * step;
		gsl_vector_set(f_step, i, temp_f_step_i);
	}

	//Compute error estimate, method from book
	for (int i = 0; i < order; i++){
		double tan0i = gsl_vector_get(tan0, i);
		double tan12i = gsl_vector_get(tan12, i);
		double temp_err_i = (tan0i - tan12i) * step / 2;
		gsl_vector_set(err, i, temp_err_i);
	}

	//Free up used memory
	gsl_vector_free(tan0);
	gsl_vector_free(tan12);
	gsl_vector_free(temp_f_val);
}




void harm_f(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    gsl_vector_set(functionDerivative, 0, gsl_vector_get(functionValue, 1));
    gsl_vector_set(functionDerivative, 1, -gsl_vector_get(functionValue, 0));
}




void rk_driver(void f(double , gsl_vector*, gsl_vector*), double a, gsl_vector* f_val_a,
	double b, gsl_vector* f_val_b, double step, double abs_acc, 
	double rel_acc, FILE* path_to_file){

	int order = f_val_a->size; //order of diff eq.
	double tol; ///tolerance
	double err;
	double norm_f;
	double pos = a;
	gsl_vector* f_step = gsl_vector_alloc(order);
	gsl_vector* current_f_val = gsl_vector_alloc(order);
	gsl_vector* f_err = gsl_vector_alloc(order);

	gsl_vector_memcpy(current_f_val,f_val_a);

	// run for full interval
	while(pos<b){
		if (path_to_file != NULL){
			fprintf(path_to_file, "%.5g \t", pos);
			for (int i = 0; i<order; i++){
				fprintf(path_to_file, "%.5g \t", gsl_vector_get(current_f_val, i));
			}
			if(f == harm_f){
			fprintf(path_to_file, "%.5g \n", sin(pos));
			}
			else{fprintf(path_to_file, "\n");}
		}

		double fin_step;
		double n_step = step;

		if(pos+n_step > b){ //overshooting interval
			n_step = b - pos;
		}

		do{
		rk_stepXY(f,pos,current_f_val,n_step,f_step,f_err);
		err = gsl_blas_dnrm2(f_err);
		norm_f = gsl_blas_dnrm2(f_step);
		tol = (norm_f * rel_acc + abs_acc) *
	                        sqrt(n_step / (b - a));
		fin_step = n_step;
		n_step *= pow(tol/err, 0.25)*0.95; // make steps smaller
		}
		while(err > tol);


		gsl_vector_memcpy(current_f_val,f_step);
		pos += fin_step;

	}

	gsl_vector_memcpy(f_val_b, f_step); //Function value at right end point

	//Free up memory
	gsl_vector_free(f_step);
 	gsl_vector_free(f_err);
	gsl_vector_free(current_f_val);
}

