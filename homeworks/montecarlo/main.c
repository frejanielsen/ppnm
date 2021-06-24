#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mc_fs.h"


int main(){


double pi = 3.1415926;
double vars[] = {0.0};
	double f1(int dim, double* vars){
		double x = vars[0];
		return exp(-x)+1;
		}

double a1[] = {0.0};
double b1[] = {5.0};
int dim = sizeof(a1)/sizeof(a1[0]);
int N = 10000;
double res = 0;
double err = 0;
printf("checking function f1 works: %f\n",f1(dim,vars));


//double res1 = plainmc(dim,f1,a1,b1,N);
//printf(" Main: plainmc called. \n Result=%f\n Actual result: 5.9933\n",res1);

plainmc(dim,f1,a1,b1,N,&res,&err);
printf("integral of exp(-x)+1 from 0 to 5 , with 10000 points:\n");
printf(" plainmc called. \n res=%f\n err=%f\n actual result: 5.9933\n",res,err);




	double f2(int dim, double* vars){
		double x = vars[0];
		double y = vars[1];
		return cos(x)*cos(y)*sin(x);
		}

double a2[] = {0.0, 0,0};
double b2[] = {pi,pi};
dim = sizeof(a2)/sizeof(a2[0]);
N = 30000;

//printf("checking function f2: %f\n",f2(dim,vars));


plainmc(dim,f2,a2,b2,N,&res,&err);
printf("integral of cos(x)*cos(y)*sin(x) from 0 to pi , with 30000 points:\n");
printf(" plainmc called. \n res=%f\n err=%f\n actual result: 0\n",res,err);



	double f3(int dim, double* vars){
		double x = vars[0];
		double y = vars[1];
		double z = vars[2];
		double pi = 3.1415926;
		return 1.0/(1-cos(x)*cos(y)*cos(z))*1.0/(pi*pi*pi);
		}

double a3[] = {0.0,0.0,0.0};
double b3[] = {pi,pi,pi};
dim = sizeof(a3)/sizeof(a3[0]);
N = 30000;
//printf("checking function f3: %f\n",f3(dim,vars));



plainmc(dim,f3,a3,b3,N,&res,&err);
printf("integral of 1/(1-cos(x)*cos(y)*cos(z)) from 0 to pi in x, y and z, with 30000 points:\n");
printf(" plainmc called. \n res=%f\n err=%f\n actual result: 1.39320393\n",res,err);











return 0;
}
