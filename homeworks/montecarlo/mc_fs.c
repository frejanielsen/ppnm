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


double corput (int n, int base){
	double q=0;
	double bk=1.0/base;
	while(n>0){
		q+=(n % base)*bk;
		n/=base;
		bk/=base;
	}
return q;
}


void multimc(int dim, double f(int dim, double* x), double* a, double* b, int N, double* res, double* err){

	double V=1;
	for (int i=0;i<dim;i++){
	        V*=b[i]-a[i];
	}

	double N_half = round(N/2);

	//double the sum variables for error estimates
	double sum1=0;
	double sum2=0;
	double x[dim];

	double fx;

	//list of primes, sets a maximum dimension on the integral, add primes if more dimensions needed.
	int primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271};

	for (int i=0;i<N_half;i++){
	        for (int j=0;j<dim;j++){
	                x[j]=a[j]+corput(i,primes[j])*(b[j]-a[j]);
	        }

	        fx = f(dim,x);
	        sum1 += fx;

		for (int j=0;j<dim;j++){
	                x[j]=a[j]+corput(i,primes[j+dim])*(b[j]-a[j]);
	        }

		fx = f(dim,x);
		sum2 += fx;
	}
	double avg=(sum2+sum1)/(2*N_half);
	*res=avg*V;
	*err=fabs(sum1-sum2)/(2*N_half)*V;

}


