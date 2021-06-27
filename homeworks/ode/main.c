#include <math.h>
#include <gsl/gsl_matrix.h>
#include "rung_kutta.h"
#include "diff_eqs.h"




int main(int argc, char * argv[]){
	printf("Results found in individual txt-files for each example run, each with corresponding plots as png files\n");

	//A
	//Solving u'' = -u
	int harmonicDimension = 2;
    gsl_vector *harmonicFunctionValueLeft = gsl_vector_alloc(harmonicDimension);
    gsl_vector *harmonicFunctionValueRight = gsl_vector_alloc(harmonicDimension);
    gsl_vector_set(harmonicFunctionValueLeft, 1, 1); //Initial value

    //Settings for ODE solver
    double leftEndpoint = 0.0;
    double rightEndPoint = 3 * M_PI;
    double absoluteAccuracy = 1e-3;
    double relativeAccuracy = 1e-3;
    double step = (rightEndPoint - leftEndpoint) / 10;

    FILE *harm_output = fopen(argv[1], "w");
    rk_driver(&harm_f, leftEndpoint, harmonicFunctionValueLeft, rightEndPoint, harmonicFunctionValueRight,
              step, absoluteAccuracy, relativeAccuracy, harm_output);
    //fprintf(harm_output, "HARMONIC OSCILATOR - PLOT IN .png FILE");
    fclose(harm_output);

    //Solving SIR-model
    int SIRDimension = 3;
    gsl_vector *SIRFunctionValueLeft = gsl_vector_alloc(SIRDimension);
    gsl_vector *SIRFunctionValueRight = gsl_vector_alloc(SIRDimension);

    leftEndpoint = 0;
    rightEndPoint = 100;

    //COVID-19 data in Denmark 12/4/2021
    double population = 5806000;
    double totalInfected = 237792;
    double recovered = 226630;
    double currentlyInfected = totalInfected - recovered;
    double vaccinated = 445566;
    double dead = 2441;
    double removed = dead + recovered + vaccinated;

    gsl_vector_set(SIRFunctionValueLeft, 0, population - currentlyInfected - removed);
    gsl_vector_set(SIRFunctionValueLeft, 1, currentlyInfected);
    gsl_vector_set(SIRFunctionValueLeft, 2, removed);

    FILE *SIR_output = fopen(argv[2], "w");
    rk_driver(&SIR_model, leftEndpoint, SIRFunctionValueLeft, rightEndPoint, SIRFunctionValueRight, step,
              absoluteAccuracy, relativeAccuracy, SIR_output);
    fclose(SIR_output);

    //SIR again, this time with doubled contact time
    gsl_vector *SIRFunctionValueLeft2 = gsl_vector_alloc(SIRDimension);
    gsl_vector *SIRFunctionValueRight2 = gsl_vector_alloc(SIRDimension);
    gsl_vector_set(SIRFunctionValueLeft2, 0, population - currentlyInfected - removed);
    gsl_vector_set(SIRFunctionValueLeft2, 1, currentlyInfected);
    gsl_vector_set(SIRFunctionValueLeft2, 2, removed);
    FILE *SIR_output2 = fopen(argv[3], "w");
    rk_driver(&SIR_model_new_contact_time, leftEndpoint, SIRFunctionValueLeft2, rightEndPoint, SIRFunctionValueRight2, step,
              absoluteAccuracy, relativeAccuracy, SIR_output2);
    fclose(SIR_output2);



return 0;
}
