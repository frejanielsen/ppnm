#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

#include "funcs.h" // include functions




int main () {

	double xmin = 0;
	double xmax = 1;

	double acc = 0.0001;
	double eps = 0.0001;

	int nrec = 0;
	int calls = 0;

	double Err = 0;



	double f(double x){
		calls++;
		return sqrt(x);
	}

	double I=integ(f,xmin,xmax,acc,eps,nrec,&Err);
	printf("integral of sqrt(x) from 0 to 1\n");
	printf("Nummerical result=%f \n Analytical result 2/3 \nFunction called %i times\n",I,calls);
	printf("Integration error %f\n",Err);
	//printf("The error estimate is somewhat pessimistic by the nature of the calculation method, this goes for the entire work.\n");

	Err = 0;
	calls = 0;
	nrec = 0;

	double f_2(double x){
		calls++;
		return 4*sqrt(1-x*x);
	}

	I=integ(f_2,xmin,xmax,acc,eps,nrec,&Err);
	printf("\n\n\nintegral of 4*sqrt(1-x²) from 0 to 1\n");
	printf("Nummerical result=%f \n Analytical result pi \nFunction called %i times\n",I,calls);
	printf("Integration error %f\n",Err);


	return 0;
}
