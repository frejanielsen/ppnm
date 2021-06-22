#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main (){


double x;
double y;
int items;
int First_Flag =0;

while(items!=EOF){
	if (First_Flag==0){
		items=scanf("%lg",&y);
		x=y;
		First_Flag=1;
	} else {
		items=scanf("%lg",&y);
		printf("x=%g, sin(x)=%g, cos(x)=%g\n",x,sin(x),cos(x));
		x=y;
	};
};


return 0;
}
