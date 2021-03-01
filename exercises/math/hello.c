#include<stdio.h>
#define PI 3.1415927

//void world(void);
#include"hello.h"

char hello[] = "file scope";
// comment line

/*
block comment
*/

double x;

void bar(void){
	char hello[] = "function scope";
	printf("hello (function scope?):%s\n",hello);
	{
		char hello[] = "block scope";
		printf("hello (block scope?):%s\n",hello);
	}
}

int main(){
	printf("hello\n");
	world();
	printf("hello (file scope?):%s\n",hello);
	bar();
	x=7;
	printf("x=%g\n",x);
return 0;
}

