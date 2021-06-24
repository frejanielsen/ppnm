#include "diff_eqs.h"
#include <math.h>
#include <gsl/gsl_vector.h>


void harm_f(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    gsl_vector_set(functionDerivative, 0, gsl_vector_get(functionValue, 1));
    gsl_vector_set(functionDerivative, 1, -gsl_vector_get(functionValue, 0));
}




void SIR_model(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    double population = 5806000; //population of denmark 2019
    double recoveryTime = 21; //Three weeks recovery
    double contactTime = 3; //Time between contacts
    double susceptible = gsl_vector_get(functionValue, 0);
    double infected = gsl_vector_get(functionValue, 1);
    double derivativeSusceptible = -(susceptible * infected) / (contactTime * population);
    double derivativeRemoved = infected / recoveryTime;
    double derivativeInfected = -derivativeSusceptible - derivativeRemoved;

    gsl_vector_set(functionDerivative, 0, derivativeSusceptible);
    gsl_vector_set(functionDerivative, 1, derivativeInfected);
    gsl_vector_set(functionDerivative, 2, derivativeRemoved);
}


void SIR_model_new_contact_time(double variable, gsl_vector *functionValue, gsl_vector *functionDerivative)
{
    double population = 5806000; //population of denmark 2019
    double recoveryTime = 21; //Three weeks recovery
    double contactTime = 6; //Time between contacts, doubled for this equation
    double susceptible = gsl_vector_get(functionValue, 0);
    double infected = gsl_vector_get(functionValue, 1);
    double derivativeSusceptible = -(susceptible * infected) / (contactTime * population);
    double derivativeRemoved = infected / recoveryTime;
    double derivativeInfected = -derivativeSusceptible - derivativeRemoved;

    gsl_vector_set(functionDerivative, 0, derivativeSusceptible);
    gsl_vector_set(functionDerivative, 1, derivativeInfected);
    gsl_vector_set(functionDerivative, 2, derivativeRemoved);
}
