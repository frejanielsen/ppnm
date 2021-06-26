#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

#include "numfs.h" // include functions

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


	/////// CC /////
/*
	nrec=0;
	calls=0;
	double I_c = integ(f_2,acc,eps,nrec);
	printf("\n\n\nIntegral of 4*sqrt(1-x²) from 0 to 1 with Clenshaw-Curtis transformation\n");
	printf("Nummerical result=%f\nAnalytical result pi \nFunction called %i times\n",I_c,calls);

	nrec=0;
	calls=0;

	double f3(double x){
		calls++;
		return 1.0/sqrt(x);
	}
	I = integ_cc(f3,acc,eps,nrec);
	printf("\n\n Integral of 1/sqrt(x) from 0 to 1 with Clenshaw-Curtis transformation:\n");
	printf("Nummerical result=%f \nAnalytical result 2 \nFunction called %i times\n",I,calls);

	nrec=0;
	calls=0;

	double f4(double x){
		calls++;
		return log(x)/sqrt(x);
	}
	I = integ_cc(f4,acc,eps,nrec);
	printf("\n\n Integral of log(x)/sqrt(x) from 0 to 1 with Clenshaw-Curtis transformation:\n");
	printf("Nummerical result=%f \nAnalytical result 2 \nFunction called %i times\n",I,calls);
*/
	return 0;
}
