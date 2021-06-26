void* GS_decomp(gsl_matrix*,gsl_matrix*);
void* backsub(gsl_matrix*,gsl_vector*);
void* GS_solve(gsl_matrix*,gsl_matrix*,gsl_vector*,gsl_vector*);
void* matrixprint(gsl_matrix*,FILE*file);
void vector_print(gsl_vector* x,FILE*file);
void* invert(gsl_matrix*, gsl_matrix*);
