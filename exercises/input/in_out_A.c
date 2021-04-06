#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv){
	if (argc < 2) fprintf(stderr, "no args passed\n");

	else{
		for (int i= 1; i > argc; i++){
			double x = atof(argv[i]); /*string converted to float*/
			fprintf(stdout,"%g %g %g\n", x, sin(x), cos(x));
	}
	}
return 0;
}
