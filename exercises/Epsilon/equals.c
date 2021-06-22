#include<stdio.h>
#include<math.h>
#include<stdlib.h>



int equal(double a, double b, double tau, double eps){
	if(abs(a-b)<tau || abs(a-b)/(abs(a)+abs(b))<eps/2){
	return 1;
	}
	else{
	return 0;
}
	}


int main(){

	int t=equal(5,6,1,2);
	printf("are they equal: %i\n",t);

	return 0;
}


