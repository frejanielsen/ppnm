#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double func_A(double x, void* params){
	double val = log(x)/(sqrt(x));
	return val;
}

double integ_A(){
	gsl_function F;
	F.function = &func_A;
	int limit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a=0, b=1, ac=1e-6, eps=1e-6, result, err;
	gsl_integration_qags(&F, a, b, ac, eps, limit, w, &result, &err);
	gsl_integration_workspace_free(w);
	return result;
}


int main(){
	FILE* out_A =fopen("out.value_exercise_A.txt","w");
	double exercise_A = integ_A();
	fprintf(out_A, "int(log(x)/sqrt(x), x=0,,1) = %g\n", exercise_A);
	printf("results found in out.value_exercise_A.txt \n");
	return 0;
}
