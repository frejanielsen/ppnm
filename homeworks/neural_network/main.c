#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>

#include "ann_fs.h"


double negative_sin(double x)
{
    return -sin(x); //Need function for derivative of cosine
}

int main(int argc, char *argv[])
{
    //PART A
    printf("A: Construct simple artificial neural network\n\n");

    int n = 5;
    int np = 20;

    //Load cosine data
    gsl_vector *xs = gsl_vector_alloc(np);
    gsl_vector *ys = gsl_vector_alloc(np);
    input_to_array(n, xs, ys, argv[1]);

    printf("Initialize neural network with %d neurons, one hidden layer \n\n", n);
    ann *network = ann_alloc(n, &cos, &negative_sin, &sin);

    printf("Training the network: \n\n");
    ann_train(network, xs, ys);

    FILE *out = fopen(argv[2], "w");

    double lower_lim = 0;
    double upper_lim = 11;

    for (double i = lower_lim; i <= upper_lim; i += 1.0 / 8)
    {
        fprintf(out, "%10g %10g %10g %10g %10g %10g %10g\n", i, ann_response(network, i), cos(i),
                ann_response_d(network, i), -sin(i),
                ann_response_i(network, 0, i), sin(i));
    }

    printf("Part A and B:\n\n");
    printf("network prediction of cos(x) and its derivative + anti-derivative can be seen in plot.png\n\n");


    fclose(out);

    ann_free(network);

    return 0;
}
