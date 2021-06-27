#include<stdio.h>
#include<math.h>

double exponential(double x){
	if(x<0)return 1/exponential(-x);
	if(x>1./8)return pow(exponential(x/2),2);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

int main(){
	for( double x = 1./16; x <= 8; x += 1./16 ){
			printf("%g %g %g\n",x,exponential(x),exp(x));
		}

return 0;
}
