#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_interp.h"


int main(){

int n=20;
double x[n], y[n];
int i=0;
while(i<n){
	x[i]=i/2.0;
	y[i]=sin(x[i]);
	i++;
}


gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear,n);
gsl_interp* cspline= gsl_interp_alloc(gsl_interp_cspline,n);
gsl_interp* polynomial= gsl_interp_alloc(gsl_interp_polynomial,n);

gsl_interp_init(linear,x,y,n);
gsl_interp_init(cspline,x,y,n);
gsl_interp_init(polynomial,x,y,n);

int z=0;
double fineness=10;


FILE* gsl_pol=fopen("out.gsl.pol.txt","w");
FILE* gsl_lin_integ=fopen("out.gsl.lin.integ.txt","w");
FILE* gsl_pol_integ=fopen("out.gsl.pol.integ.txt","w");

for (z=0;z<fineness*x[n-1];z++){
	printf("%10g %10g \n",z/fineness,gsl_interp_eval(linear,x,y,z/fineness,NULL));
	fprintf(gsl_pol,"%10g %10g \n",z/fineness,gsl_interp_eval(polynomial,x,y,z/fineness,NULL));
	fprintf(gsl_lin_integ,"%10g %10g \n",z/fineness,gsl_interp_eval_integ(linear,x,y,x[0],z/fineness,NULL));
	fprintf(gsl_pol_integ,"%10g %10g \n",z/fineness,gsl_interp_eval_integ(polynomial,x,y,x[0],z/fineness,NULL));
}






return 0;
}
