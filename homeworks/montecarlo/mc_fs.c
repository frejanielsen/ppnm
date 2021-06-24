#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

/*
complex double plainmc(int dim, double f(int dim, double* x), double* a, double* b,
	int N){

        double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0,sum2=0,x[dim];
        for(int i=0;i<N;i++){
                for(int i=0;i<dim;i++)x[i]=a[i]+RAND_MAX*(b[i]-a[i]);
                double fx=f(dim,x); sum+=fx; sum2+=fx*fx;
                }
        double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
        complex result=mean*V+I*sigma*V/sqrt(N);
        return result;
}
*/

void plainmc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* err){


	double V=1;
	for (int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}

	double sum=0;
	double sum2=0;
	double x[dim];
	unsigned int rand_seed=3;
	double fx;

	for (int i=0;i<N;i++){
		for (int j=0;j<dim;j++){
			x[j]=a[j]+rand_r(&rand_seed)*1.0/RAND_MAX*(b[j]-a[j]);
		}

		fx=f(dim,x);
		sum+=fx;
		sum2+=fx*fx;
	}
	double avg=sum/N;
	*res=avg*V;
	*err=sqrt(sum2/N-avg*avg)*V/sqrt(N);

}
