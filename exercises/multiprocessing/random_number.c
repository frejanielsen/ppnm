#include <stdlib.h>
#include "random_number.h"

double random_number(unsigned int *seed)
{
    double max_random = (double) RAND_MAX;
    double random_number = (double) rand_r(seed);

    return 2*random_number/max_random-1; //rescale between -1 and 1
}
