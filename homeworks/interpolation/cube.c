#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>


typedef struct {int n; double* x, *y, *b, *c, *d;} cubicspline;
	cubicspline* cubicspline_alloc(int n, double* x, double* y){
	cubicspline *s = (cubicspline*)malloc(sizeof(cubicspline));
	s->b = (double*)malloc((n)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;

	int i=0;
	for(i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	double dy[n-1], dx[n-1];
	for (i=0;i<n-1;i++){
		dx[i]=x[i+1]-x[i];
		assert(dx[i]>0);
		dy[i]=(y[i+1]-y[i])/dx[i];
	}
	double D[n], Q[n-1], B[n];
	D[0]=2;
	Q[0]=1;
	for(i=0;i<n-2;i++){
		D[i+1]=2*dx[i]/dx[i+1]+2;
		Q[i+1]=dx[i]/dx[i+1];
		B[i+1]=3*(dy[i]+dy[i+1]*dx[i]/dx[i+1]);
	}
	D[n-1]=2;
	B[0]=3*dy[0];
	B[n-1]=3*dy[n-2];
	for(i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];
	}
	s->b[n-1]=B[n-1]/D[n-1];
	for(i=n-2;i>=0;i--){
		s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	}
	for(i=0;i<n-1;i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*dy[i])/dx[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*dy[i])/dx[i]/dx[i];
	}
	return s;
}



int binsearch(int n, double* x, double input){
assert(n>1 && input<=x[n-1] && input>=x[0]);
int i=0, j=n-1, m=0; // setting up for binary search for intervals
while(j-i>1){
	m=(i+j)/2;

	if(input<x[m]){
		j=m;
	}else{
		i=m;
	}
}
return i;
}


double cubicspline_eval(cubicspline *s, double input){
	int i = binsearch(s->n,s->x,input);
	double dx=input-s->x[i];
	double output=s->y[i]+dx*(s->b[i]+dx*(s->c[i]+dx*s->d[i]));
return output;
}


double cubicspline_integ(cubicspline *s, double input){
	int i = binsearch(s->n,s->x,input);
	double dx=0;
	double output=0;
	for (int p=1;p<=i;p++){
		dx=s->x[p]-s->x[p-1];
		output+=dx*(s->y[p-1]+dx*(s->b[p-1]/2+dx*(s->c[p-1]/3+s->d[p-1]/4*dx)));
	}
	dx=input-s->x[i];
	output+=dx*(s->y[i]+dx*(s->b[i]/2+dx*(s->c[i]/3+s->d[i]/4*dx)));

return output;
}



void cubicspline_free(cubicspline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}





int main(){

int n=20;

double x[n],y[n];

FILE* my_out_stream=fopen("out.xydata.txt","w");

for(int i=0;i<n;i++){
x[i]=i/2.0;
y[i]=sin(x[i]);
//y[i]=i*i;
fprintf(my_out_stream,"%10g %10g \n",x[i],y[i]);
}

cubicspline* s = cubicspline_alloc(n,x,y);


FILE* out_cubic_integ = fopen("out.cubic.integ.txt","w");


int z=0;
double fineness=10;

while(z<=fineness*s->x[n-1]){
	printf("%10g	%10g\n",z/fineness,cubicspline_eval(s,z/fineness));
	fprintf(out_cubic_integ,"%10g %10g\n",z/fineness,cubicspline_integ(s,z/fineness));
	z++;
}


cubicspline_free(s);

return 0;
}
