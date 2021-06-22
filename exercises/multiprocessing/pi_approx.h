#ifndef UNI_APPROXIMATEPI_H
#define UNI_APPROXIMATEPI_H

#include <pthread.h>

void *place_points(void *place_points_struct_input);

void approximate_pi_multi(double number_of_points);

#endif //UNI_APPROXIMATEPI_H
