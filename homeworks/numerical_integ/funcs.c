#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "funcs.h"

double integ_reuse (double f(double),double xmin, double xmax,
	double acc, double eps, double f2, double f3,  int nrec){
	assert(nrec<50000);

	double f1 = f(xmin + (xmax - xmin) / 6); //function evaluations not reused
	double f4 = f(xmin + 5 * (xmax - xmin) / 6);
	double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * (xmax - xmin); //result of 4. order
	double q = (f1 + f2 + f3 + f4) / 4 * (xmax - xmin);	//result of 2. order
	double tol = acc+eps*fabs(Q);
	double err = fabs(Q - q);

	double* Err = 0;
	if(err < tol){
		*Err = sqrt((*Err) * (*Err) + err * err);
		return Q;
	} else {
		nrec++;
		double Qi = integ_reuse(f,xmin,(xmax + xmin)/2,acc/sqrt(2.0),eps,f1,f2,nrec);
		double Qii = integ_reuse(f,(xmax + xmin)/2,xmax,acc/sqrt(2.0),eps,f3,f4,nrec);
		return Qi+Qii;
	}
}

double integ_cc(double f(double), double acc, double eps, int nrec){
	double xmin = 0;
	double xmax = 3.14159;
	double F(double x){
		return F(cos(x))*sin(x);
	}
	double f2 = F(xmin+(xmax-xmin)/3);
	double f3 = F(xmin+2*(xmax-xmin)/3);
	return integ_reuse(F,xmin,xmax,acc,eps,f2,f3,nrec);
}

double integ(double f(double),double xmin, double xmax, double acc, double eps, int nrec){

	double f2 = 0;
	double f3 = 0;
	double xMax = xmax;
	double xMin = xmin;

	double result=0;

	if (isinf(xmax)==1 && isinf(xmin)==-1){ // in infinate
		xMax = 1;
		xMin = -1;

		double F (double x){
			return f(x/(1.0 - x * x)) * (1.0 + x * x)/((1.0 - x * x) * (1.0 - x * x));
		}

		f2 = F(xMin + (xMax - xMin)/3);
		f3 = F(xMin + 2*(xMax - xMin)/3);

		result = integ_reuse(F,xMin,xMax,acc,eps,f2,f3,nrec);
	}

	if (isinf(xmax)==1 && isinf(xmin)==0) {
		xMin = 0;
		xMax = 1;
		double F (double x){
	                return f(xmin + x/(1.0 - x))/((1.0 - x)*(1.0 - x));
	        }

		f2 = F(xMin+(xMax-xMin)/3);
	        f3 = F(xMin+2*(xMax-xMin)/3);


	        result = integ_reuse(F,xMin,xMax,acc,eps,f2,f3,nrec);
	}


	if (isinf(xmin)==-1 && isinf(xmax)==0){
		xMin = 0;
		xMax = 1;

		double F (double x){
	                return f(xmax-(1-x)/x)/(x*x);
	        }

		f2 = F(xMin+(xMax-xMin)/3);
	        f3 = F(xMin+2*(xMax-xMin)/3);


	        result = integ_reuse(F,xMin,xMax,acc,eps,f2,f3,nrec);
	}

	if (isinf(xmin)==0 && isinf(xmax)==0){
	//        i=3;
	        xMin=xmin;
	        xMax=xmax;

	        double F (double x){
	                return f(x);
	        }

		f2=F(xMin+(xMax-xMin)/3);
	        f3=F(xMin+2*(xMax-xMin)/3);


	        result = integ_reuse(F,xMin,xMax,acc,eps,f2,f3,nrec);
	}


	return result;
}



