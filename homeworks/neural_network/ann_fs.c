#include <time.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <float.h>
#include "ann_fs.h"

//// lin alg functions /////
double random_number(unsigned int *seed)
{
    double maxRand = (double) RAND_MAX;           // Maximum random number, cast to double
    double randNum = (double) rand_r(seed);     // Generate pseudo-random number from seed, cast to double
    return randNum / maxRand;
}

void printvector(char *s, gsl_vector *v)
{
    printf("%s\n", s);
    for (int i = 0; i < v->size; i++)
    {
        printf("%10g ", gsl_vector_get(v, i));
    }
    printf("\n");
}

void print_matrix(int n, gsl_matrix *A, char *s)
{
    printf("\n%s\n", s);
    for (int r = 0; r < n; r++)
    {
        gsl_vector_view A_r = gsl_matrix_row(A,r);
        gsl_vector *v = &A_r.vector;
        for (int i = 0; i < v->size; i++)
        {
            if (gsl_vector_get(v,i) > 1e-10)
            {
                printf("%10g\t", gsl_vector_get(v,i));
            }
            else
            { printf("%10g\t", 0.0); }
        }
        printf("\n");
    }
}

void symmetric(gsl_matrix *A, unsigned int *seed)
{
    for (int r = 0; r < (A->size1); r++)
    {
        gsl_matrix_set(A,r,r, random_number(seed));

        for (int c = r + 1; c < A->size2; c++)
        {
            double e = random_number(seed);
            gsl_matrix_set(A,r,c,e);
            gsl_matrix_set(A,c,r,e);
        }
    }

}
////// Minimization functions //////
#define MY_DOUBLE_EPSILON 2.22045e-10


void gradient(double f(gsl_vector *), gsl_vector *m, gsl_vector *g)
{
    long double step_s = MY_DOUBLE_EPSILON; //Long double?
    double f_val = f(m);
    int dim = m->size;

    for (int i = 0; i < dim; i++)
    {
        double step;
        double m_i = gsl_vector_get(m, i);

        if (fabs(m_i) < step_s)
        {
            step = step_s;
        }
        else
        {
            step = fabs(m_i) * step_s;
        }

        gsl_vector_set(m, i, m_i + step);
        gsl_vector_set(g, i, (f(m) - f_val) / step);
        gsl_vector_set(m, i, m_i - step);
    }
}

void multimc(double f(gsl_vector *), gsl_vector *m, double tol)
{
    double step_s = MY_DOUBLE_EPSILON;
    int dim = m->size;
    int nst = 0;
    int nsc = 0;
    int nr = 0;

    gsl_matrix *H = gsl_matrix_alloc(dim,dim);
    gsl_matrix *I = gsl_matrix_alloc(dim,dim);
    gsl_matrix_set_identity(H);
    gsl_matrix_set_identity(I);

    //Allocate memory for needed parts
    gsl_vector *g_val = gsl_vector_alloc(dim);
    gsl_vector *n_g_val = gsl_vector_alloc(dim);
    gsl_vector *news = gsl_vector_alloc(dim);
    gsl_vector *n_m = gsl_vector_alloc(dim);
    gsl_vector *sol = gsl_vector_alloc(dim);
    gsl_vector *sol_c = gsl_vector_alloc(dim);
    gsl_vector *b = gsl_vector_alloc(dim);

    gradient(f, m, g_val);
    double f_val = f(m);
    double n_f_val;

    while (nst < 1e4)
    {
        nst++;
        gsl_blas_dgemv(CblasNoTrans, -1, H, g_val, 0, news);
        if (gsl_blas_dnrm2(news) < step_s * gsl_blas_dnrm2(m))
        {
            fprintf(stderr, "Quasi_newton_method: |dx| < stepSize*|x|\n");
            break;
        }
        if (gsl_blas_dnrm2(g_val) < tol)
        {
            fprintf(stderr, "Quasi_newton_method: |grad| < accuracy\n");
            break;
        }

        double scale = 1;

        while (1) //while(1) means running until explict break
        {
            gsl_vector_memcpy(n_m, m);
            gsl_vector_add(n_m, news);
            n_f_val = f(n_m);
            double s_g;
            gsl_blas_ddot(news, g_val, &s_g);

            if (n_f_val < f_val + 0.01 * s_g)
            {
                nsc++;
                break;
            }
            if (scale < step_s)
            {
                nr++;
                gsl_matrix_set_identity(H);
                break;
            }
            scale *= 0.5;
            gsl_vector_scale(news, 0.5);
        }

        gradient(f, n_m, n_g_val);
        gsl_vector_memcpy(sol, n_g_val);
        gsl_blas_daxpy(-1, g_val, sol);
        gsl_vector_memcpy(sol_c, news);
        gsl_blas_dgemv(CblasNoTrans, -1, H, sol, 1, sol_c);

        gsl_matrix *sc_sc = gsl_matrix_calloc(dim, dim); //u*u^t
        gsl_blas_dsyr(CblasUpper, 1.0, sol_c, sc_sc);
        double sct; //u^T*y
        gsl_blas_ddot(sol_c, sol, &sc_sc);
        if (fabs(sct) > 1e-12)
        {
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0 / sct, sc_sc, I, 1.0, H);
        }

        gsl_vector_memcpy(m, n_m);
        gsl_vector_memcpy(g_val, n_g_val);
        f_val = n_f_val;
    }

    //Free allocated memory
    gsl_matrix_free(H);
    gsl_matrix_free(I);
    gsl_vector_free(g_val);
    gsl_vector_free(n_g_val);
    gsl_vector_free(news);
    gsl_vector_free(n_m);
    gsl_vector_free(sol);
    gsl_vector_free(sol_c);
    gsl_vector_free(b);

    fprintf(stderr,
            "Quasi_newton_method: \n amount of steps = %i \n amount of scales = %i, \n amount of matrix resets = %i \n  f(x) = %.1e\n\n",
            nst, nsc, nr, f_val);
}




////// Neural network functions /////


ann *ann_alloc(int n, double (*f)(double),double(*df)(double),double(*fi)(double))
{
    int np = 3;
    ann *nn = (ann *) malloc(sizeof(ann));
    nn->params = gsl_vector_alloc(n * np);
    nn->f = f;
    nn->df = df;
    nn->fi = fi;
    nn->n = n;

    return nn;
}

double ann_response(ann *network, double ep)
{
    int np = 2;
    int n = network->n;
    double res = 0;

    for (int i = 0; i < n; ++i)
    {
        double nsh = gsl_vector_get(network->params, np * i);
        double nsc = gsl_vector_get(network->params, np * i + 1);
        double ned = gsl_vector_get(network->params, np * i + 2);

        res += (network->f((ep - nsh) / nsc)) * ned;
    }
    return res;
}
double ann_response_d(ann *network, double ep)
{
    int np = 2;
    int n = network->n;
    double res = 0;

    for (int i = 0; i < n; ++i)
    {
        double nsh = gsl_vector_get(network->params, np * i);
        double nsc = gsl_vector_get(network->params, np * i + 1);
        double ned = gsl_vector_get(network->params, np * i + 2);

        res += (network->df((ep - nsh) / nsc)) * ned;
    }
    return res;
}
double ann_response_i(ann *network, double rp, double lp)
{
    int np = 2;
    int n = network->n;
    double res = 0;

    for (int i = 0; i < n; ++i)
    {
        double nsh = gsl_vector_get(network->params, np * i);
        double nsc = gsl_vector_get(network->params, np * i + 1);
        double ned = gsl_vector_get(network->params, np * i + 2);

        res += (((network->fi((lp - nsh) / nsc)) * ned)-((network->fi((rp - nsh) / nsc)) * ned));
    }
    return res;
}

void ann_train(ann *network, gsl_vector *data, gsl_vector *labels)
{
    unsigned int seed = time(NULL);
    int np = 3;
    int n_p = data->size;
    int n = network->n;

    double cost(gsl_vector *nextp)
    {
        int n = network->n;
        gsl_vector *updatedp = gsl_vector_alloc(np * n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < np; ++j)
            {
                gsl_vector_set(updatedp, np * i + j,
                               gsl_vector_get(nextp, np * i + j));
            }
        }

        network->params = updatedp; //Update parameters
        double cos = 0;

        for (int i = 0; i < n_p; ++i)
        {
            double ep = gsl_vector_get(data, i);
            double p_label = gsl_vector_get(labels, i);
            double res = ann_response(network, ep);
            cos += (res - p_label) * (res - p_label);
        }

        return cos;
    }

    double tol = 1e-5;
    gsl_vector *lp = gsl_vector_alloc(n * np);

    double a = -5;
    double b = 5;

    for (int i = 0; i < n; i++)
    {
        double weight = 1.000001;
        double scale = 1;
        double shift = -5+(b-a)*i/(n-1);

        gsl_vector_set(lp, 3*i,1);
        gsl_vector_set(lp, 3*i+1,scale);
        gsl_vector_set(lp, 3*i+2,weight);
    }

    multimc(cost, lp, tol);

    network->params = lp;

}


void ann_free(ann *network)
{
    gsl_vector_free(network->params);
    free(network);
}
