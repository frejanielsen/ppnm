#include <unistd.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>

#include "linfs.h"
#include "rootfs.h"
#include "odefs.h"




void f1(gsl_vector* x, gsl_vector* fx){
	gsl_blas_dcopy(x,fx);
}
void f2(gsl_vector* x, gsl_vector* fx){
        double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);

        gsl_vector_set(fx,0,sin(x0)*cos(x1));
        gsl_vector_set(fx,1,cos(x0)*sin(x0));
}
void fRv(gsl_vector* x, gsl_vector* fx){
        double x0=gsl_vector_get(x,0);
        double x1=gsl_vector_get(x,1);

        gsl_vector_set(fx,0,2*(1-x0)*(-1)+200*(x1-x0*x0)*(-2*x0));
        gsl_vector_set(fx,1,200*(x1-x0*x0));
}



int main(){

unsigned int rand_seed=3;
FILE* out = fopen("out_root.txt","w");
fprintf(out,"A: \n\n");

////// 1D function /////
int dim = 1;
fprintf(out,"1d function:\n");

gsl_vector* x_1 = gsl_vector_alloc(dim);
gsl_vector* fx_1 = gsl_vector_alloc(dim);

for (int i=0;i<dim;i++){
	gsl_vector_set(x_1,i,rand_r(&rand_seed)*1.0/RAND_MAX*15);
}

fprintf(out,"start x:\n");
vector_print(x_1,out);
fprintf(out,"f(X): \n");
f1(x_1,fx_1);
vector_print(fx_1,out);
double eps = 0.0001;
fprintf(out,"root found from algorithm:\n");
newton(f1,x_1,eps);
vector_print(x_1,out);
fprintf(out,"f(x_root):\n");
vector_print(fx_1,out);
fprintf(out,"with accuracy: %10f\n\n\n",eps);

///// 2D function /////
dim = 2;
fprintf(out,"2d function:\n");

gsl_vector* x_2 = gsl_vector_alloc(dim);
gsl_vector* fx_2 = gsl_vector_alloc(dim);

for (int i=0;i<dim;i++){
	gsl_vector_set(x_2,i,rand_r(&rand_seed)*1.0/RAND_MAX*15);
}

fprintf(out,"start x:\n");
vector_print(x_2,out);
fprintf(out,"f(X): \n");
f2(x_2,fx_2);
vector_print(fx_2,out);
eps = 0.0001;
fprintf(out,"root found from algorithm:\n");
newton(f2,x_2,eps);
vector_print(x_2,out);
fprintf(out,"f(x_root):\n");
vector_print(fx_2,out);
fprintf(out,"with accuracy: %10f\n\n\n",eps);


///// RV function /////
dim = 2;
fprintf(out,"the gradient of the Rosenbrock's valley function:\n");

gsl_vector* x_3 = gsl_vector_alloc(dim);
gsl_vector* fx_3 = gsl_vector_alloc(dim);

for (int i=0;i<dim;i++){
	gsl_vector_set(x_3,i,rand_r(&rand_seed)*1.0/RAND_MAX*15);
}

fprintf(out,"start x:\n");
vector_print(x_3,out);
fprintf(out,"f(X): \n");
fRv(x_3,fx_3);
vector_print(fx_3,out);
eps = 0.00005;

fprintf(out,"root found from algorithm:\n x:\n");
newton(fRv,x_3,eps);
vector_print(x_3,out);
fprintf(out,"f(x_root):\n");
vector_print(fx_3,out);
fprintf(out,"with accuracy: %10f\n\n\n",eps);


/////  Hydrogen  //////






return 0;
}
