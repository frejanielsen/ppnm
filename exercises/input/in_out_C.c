#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main (int argc, char** argv){


double x;
double y;
int items;
int First_Flag =0;

FILE* in=fopen(argv[1],"r");
FILE* Out=fopen(argv[2],"w");

while(items!=EOF){
	if (First_Flag==0){
		items=fscanf(in,"%lg",&y);
		x=y;
		First_Flag=1;
	} else {
		items=fscanf(in,"%lg",&y);
		fprintf(Out,"x=%g, sin(x)=%g, cos(x)=%g\n",x,sin(x),cos(x));
		x=y;
	};
};

fclose(in);
fclose(Out);
return 0;
}
