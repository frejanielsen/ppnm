#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
	int n=10;
	double x,y=1;
	complex z;
	complex im=csqrt(-1);
	complex zz=csqrt(im);
	z=csqrt(-2);
	printf("sqrt(-2)=%g +I%g\n",creal(z),cimag(z));
	printf("sqrt(im)=%g +I%g\n",creal(zz),cimag(zz));
	
return 0;
}

