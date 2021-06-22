#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int binsearch(int n, double * x, double * y, double z){/*from interpolationchapter*/
	assert(n>1 && z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1; /*binary search:*/
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}


double linterp(int n, double * x, double * y, double z){
	int i =  binsearch(n,x,y,z); /*use binsearch*/
	double dy = y[i+1]-y[i], dx=x[i+1]-x[i]; assert(dx>0);
	return y[i]+dy/dx*(z-x[i]);
}


double linterp_integ(int n, double* x, double* y, double z){
	int i = binsearch(n,x,y,z);
	double dx=x[i+1]-x[i]; assert(dx>0);
	int p_i = 1;
	double output = 0;
	while (p_i <= i){
		output += y[p_i-1] * (x[p_i]-x[p_i-1]) + 0.5 * (y[p_i]-y[p_i-1])*(x[p_i]-x[p_i-1]); /*intg over x[p_i-1] and x[p_i]*/
		p_i++;
	}
	double y_z = linterp(n,x,y,z);
	output += y[i] * (z-x[i]) + 0.5 * (y_z-y[i]) * (z-x[i]);
	return output;
}



int main(){
	int n = 100;
	//save output in file...
	double x[n], y[n], Y[n];

	FILE* my_output = fopen("out.xydata.txt","w");


	for (int i = 0; i < n; i++){
		x[i] = i/2.0;
		y[i] = sin(x[i]);
		Y[i] = cos(x[i])-1;
		fprintf(my_output,"%10g %10g \n", x[i],y[i]);
	}

	FILE* my_output_points = fopen("out.integ.exact.txt","w");
	for(int i = 0; i < n; i++){
		fprintf(my_output_points,"%10g	%10g\n",x[i],Y[i]);
	}

	int z = 0;
	double f = 10.0;
	while(z < f*x[n-1]){
		printf("%10g %10g\n",z/f,linterp(n,x,y,z/f));
		z++;
	}


	FILE* my_output_integ = fopen("out.xyinteg.txt","w");
	z = 0;

	while(z<=f*x[n-1]){
		fprintf(my_output_integ,"%10g	%10g\n",z/f,linterp_integ(n,x,y,z/f));
		z++;
	}

	double xs[n/2], ys[n/2];
	for (int i=0;i<n/2;i++){
		ys[i]=y[2*i];
		xs[i]=x[2*i];
	}


	FILE* my_output_half = fopen("out.xyinteg.half.txt","w");

	z=0;

	while(z<=f*xs[n/2-1]){
		fprintf(my_output_half,"%10g	%10g\n",z/f,linterp_integ(n/2,xs,ys,z/f));
		z++;
	}




return 0;
}

