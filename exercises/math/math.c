#include<math.h>
#include<complex.h>
#include<stdio.h>

# define M_E  2.7182818284590452354 /* e */
# define M_PI 3.1415926535897932384 /* pi */




int main(){

	// 1)
	// use functions, all from https://en.wikipedia.org/wiki/C_mathematical_functions#Type-generic_functions
	double complex gamma = tgamma(5);
	double complex bessel = j1(0.5);
	double complex sqrtm2 = csqrt(-2);
	double complex ePowI = cexp(I);
	double complex ePowIPi = cexp(I * M_PI);
	double complex iPowE = cpow(I, M_E);
	double complex iPowI = cpow(I, I);

	//use "crealf" to return real part of the number.
	printf("famma(5): %g + %gi \n", crealf(gamma)  , cimagf(gamma)  );
	printf("j1(5): %g + %gi \n", crealf(bessel) , cimagf(bessel) );
	printf("complex-sqrt(-2): %g + %gi \n", crealf(sqrtm2) , cimagf(sqrtm2) );
	printf("e^{i}: %g + %gi \n", crealf(ePowI)  , cimagf(ePowI)  );
	printf("e^{i*pi}: %g + %gi \n", crealf(ePowIPi), cimagf(ePowIPi));
	printf("i^{e}: %g + %gi \n", crealf(iPowE)  , cimagf(iPowE)  );
	printf("i^{i}: %g + %gi \n", crealf(iPowI)  , cimagf(iPowE)  ); 


	// 2)
	float xf = 1.f/9 ;
	double xd = 1./9  ;
	long double xld = 1.L/9 ;

	printf("float: %.25g  \n", xf  );
	printf("double: %.25lg \n", xd  );
	printf("long double: %.25Lg \n", xld );


return 0;
}

