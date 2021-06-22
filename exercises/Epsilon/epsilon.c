#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
/*1.i*/
	int i = 300000;	/*maximum value of i*/

	while(i+1>i) {i++;}
	printf("my maximum int: %i\n",i);

	int n = 1;
	for(int i=1; i<i+1; i++){n=i;}
	printf("my max int using for loop: %i\n", n+1);

	do{i++;} while(i+1>i);
	printf("my max int using do while loop: %i\n",i);

/*1.ii*/
	while(i>i-1){i--;}
	printf("Min int = %i\n",i);

	n = 1;
	for (int i=1;i>i-1;i--){n = i;}
	printf("min int usinf for loop: %i\n",n-1);

	do{i--;} while(i<i-1);
	printf("min int using do while loop: %i\n",-(i-1));

/*1.iii*/
	double x = 1;
	long double xl = 1;

	while(1+x!=1){x/=2;}
	x*=2;
	printf("min abs double: %g\n",x);

	float f = 1;
	for (float xf=1; 1+xf!=1; xf/=2){f=xf;}
	f*=2;
	printf("min abs float: %g\n", f);

	do{xl/=2;} while(1+xl!=1);
	xl*=2;
	printf("min abs long double: %Lg\n",xl);



/*2*/
	int max = INT_MAX/2;
	int c = 1;
	float g = 0;

	while(c<max+1){
	g += 1.0f/c;
	c++;
	}

	printf("sum up: %f\n",g);
	g=0;

	while(c>0){
	g += 1.0f/c;
	c--;
	}

	printf("sum down: %f\n",g);
	g=0;

	double g_d = 0;

	while(c<max+1){
	g_d += 1.0/c;
	c++;
	}

	printf("sum up: %g\n",g_d);
	g_d=0;

	while(c>0){
	g_d += 1.0/c;
	c--;
	}


	printf("sum down: %g\n",g_d);


return 0;
}

