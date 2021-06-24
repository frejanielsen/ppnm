#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "min_fs.h"



int main(){
/*
	int dim1=1;
	gsl_vector* x1=gsl_vector_alloc(dim1);
	gsl_vector_set(x1,0,3);
	double eps1=0.5;

	double f1(gsl_vector* x){
		double X=gsl_vector_get(x,0);
		return X*X;
	}

	printf("minimize x^2\n");
	printf("starting with x= %f\n", gsl_vector_get(x1,0));
	for (int i =1;i<dim1;i++){
		printf("           %10f\n",gsl_vector_get(x1,i));}

	printf("Gives: f1(x)=%10f\n",f1(x1));

	qnewton(f1,x1,eps1);
//	gsl_vector* dF=gsl_vector_alloc(dim1);
//	gradient(f1,x1,dF);
	printf("minimization done:\n x=%10f\n",gsl_vector_get(x1,0));
	for (int i =1;i<dim1;i++){
			printf("   %10f\n",gsl_vector_get(x1,i));}
	printf("gives: f(x)=%10f",f1(x1));
	printf("\n\n\n");
*/
/////////// A //////////////

	int dim2=2;
	gsl_vector* x2 = gsl_vector_alloc(dim2);
	gsl_vector_set(x2,0,10);
	gsl_vector_set(x2,1,10);
	double eps2=0.01;

	double f2(gsl_vector* x){
		double X = gsl_vector_get(x,0);
		double Y = gsl_vector_get(x,1);
		return (1-X)*(1-X)+100*(Y-X*X)*(Y-X*X);
	}

	printf("minimize (1-X)*(1-X)+100*(Y-X*X)*(Y-X*X)\n");
	printf("starting with x= %f\n", gsl_vector_get(x2,0));
	for (int i =1;i<dim2;i++){
		printf("           %10f\n",gsl_vector_get(x2,i));}

	printf("Gives: f2(x,y)=%10f\n",f2(x2));

	qnewton(f2,x2,eps2);

	printf("minimization done:\n x=%10f\n",gsl_vector_get(x2,0));
	for (int i =1;i<dim2;i++){
			printf("   %10f\n",gsl_vector_get(x2,i));}
	printf("gives: f2(x,y)=%10f",f2(x2));
	printf("\n\n\n");


	gsl_vector_set(x2,0,3);
	gsl_vector_set(x2,1,3);
	double eps3=0.05;

	double f3(gsl_vector* x){
		double X = gsl_vector_get(x,0);
		double Y = gsl_vector_get(x,1);
		return (X*X+Y-11.0)*(X*X+Y-11.0)+(X+Y*Y-7.0)*(X+Y*Y-7.0);
	}

	printf("minimize (X*X+Y-11.0)*(X*X+Y-11.0)+(X+Y*Y-7.0)*(X+Y*Y-7.0)\n");
	printf("starting with x= %f\n", gsl_vector_get(x2,0));
	for (int i =1;i<dim2;i++){
		printf("           %10f\n",gsl_vector_get(x2,i));}

	printf("Gives: f3(x,y)=%10f\n",f3(x2));

	qnewton(f3,x2,eps3);

	printf("minimization done:\n x=%10f\n",gsl_vector_get(x2,0));
	for (int i =1;i<dim2;i++){
			printf("   %10f\n",gsl_vector_get(x2,i));}
	printf("gives: f3(x,y)=%10f",f3(x2));
	printf("\n\n\n");


///////////// B /////////////
double breit_wigner(double E,gsl_vector* X){
                double A=gsl_vector_get(X,0);
                double m=gsl_vector_get(X,1);
                double Gamma=gsl_vector_get(X,2);

                return A/((E-m)*(E-m)+Gamma*Gamma/4.0);
        }


double f4(gsl_vector* X){  // mimimize this.
	int ndata=20;
	double Es[]={101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159};
	double c_sec[]={-0.25, -0.3, -0.15, -1.71, 0.81, 0.65, -0.91, 0.91, 0.96, -2.52, -1.01, 2.01, 4.83, 4.58, 1.26, 0.45, 0.15, -0.9,1 -0.81, -1.41, 1.36, 0.5, -0.45, 1.61, -2.21, -1.86, 1.76, -0.5};
	double errs[]={2, 2, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.9, 0.9, 0.9};

	double s=0;
	double bw;

	for (int i=0;i<ndata;i++){
	bw=breit_wigner(Es[i],X);
	double c=c_sec[i];
	double err=errs[i];
	s+=(bw-c)*(bw-c)/(err*err);
	}
	return s;
}

	int dim4 = 3;
	gsl_vector* x4 = gsl_vector_alloc(dim4);
	gsl_vector_set(x4,0,1e3);
	gsl_vector_set(x4,1,1);
	gsl_vector_set(x4,2,1e-1);
	double eps4 = 1e-10;

	printf("minimize A/((E-m)*(E-m)+Gamma*Gamma/4.0)\n");
	printf("starting with x= %f\n", gsl_vector_get(x4,0));
	for (int i =1;i<dim4;i++){
		printf("           %10f\n",gsl_vector_get(x4,i));}

	printf("Gives: f4(E,gamma,A)=%10f\n",f4(x4));

	qnewton(f4,x4,eps4);

	printf("minimization done:\n x=%10f\n",gsl_vector_get(x4,0));
	for (int i =1;i<dim4;i++){
			printf("   %10f\n",gsl_vector_get(x4,i));}
	printf("gives: f4(E,gamma,A)=%10f",f4(x4));
	printf("\n\n\n");

	FILE* stream = fopen("bw_fit.txt","w");
	double bw;
	for(int i=0;i<60*10;i++){
		bw=breit_wigner(100+i*1.0/10,x4);
		fprintf(stream,"%10f  %10f\n",100+i*1.0/10,bw);
	}






return 0;
}
