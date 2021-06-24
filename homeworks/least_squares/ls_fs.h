void* GS_decomp(gsl_matrix*,gsl_matrix*);
void* backsub(gsl_matrix*,gsl_vector*);
void* GS_solve(gsl_matrix*,gsl_matrix*,gsl_vector*,gsl_vector*);
void* matrixprint(gsl_matrix*);
void* fit(gsl_matrix*, double(int,double), int,gsl_vector*, gsl_matrix*);
void* invert(gsl_matrix*, gsl_matrix*);
double func(int,double);
