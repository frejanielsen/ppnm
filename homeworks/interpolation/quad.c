#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>


typedef struct {int n; double* x, *y, *b, *c;} quadspline;
	quadspline* quadspline_alloc(int n, double* x, double* y){
	quadspline *s = (quadspline*)malloc(sizeof(quadspline));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
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
		dy[i]=(y[i+1]-y[i])/dx[i];
	}
	s->c[0]=0;
	for (i=0;i<n-2;i++){
		s->c[i+1]=(dy[i+1]-dy[i]-s->c[i]*dx[i])/dx[i+1];
	}
	s->c[n-2]/=2;
	for (i=n-3;i>=0;i--){
		s->c[i]=(dy[i+1]-dy[i]-s->c[i+1]*dx[i+1])/dx[i];
	}
	for (i=0;i<n-1;i++){
		s->b[i]=dy[i]-s->c[i]*dx[i];

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


double quadspline_eval(quadspline *s, double input){
	int i = binsearch(s->n,s->x,input);
	double dx=input-s->x[i];
	double output=s->y[i]+dx*(s->b[i]+dx*s->c[i]);
return output;
}


double quadspline_integ(quadspline *s, double input){
	int i = binsearch(s->n,s->x,input);
	double output = 0;
	double dx=0;
	for (int p=1;p<=i;p++){
		dx=s->x[p]-s->x[p-1];
		output+=dx*(s->y[p-1]+dx*(s->b[p-1]/2+dx*s->c[p-1]/3));
	}
	dx=input-s->x[i];
	output+=dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));
	

return output;
}


void quadspline_free(quadspline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
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

quadspline* s = quadspline_alloc(n,x,y);

FILE* out_quad_integ = fopen("out.quad.integ.txt","w");

int z=0;
double fineness=10;

while(z<=fineness*s->x[n-1]){
printf("%10g	%10g\n",z/fineness,quadspline_eval(s,z/fineness));
fprintf(out_quad_integ,"%10g %10g\n",z/fineness,quadspline_integ(s,z/fineness));
z++;
}





quadspline_free(s);

return 0;
}
