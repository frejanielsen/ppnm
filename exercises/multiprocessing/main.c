#include <pthread.h>
#include <stdio.h>
#include "pi_approx.h"

int main()
{
    double number_of_points = 1e7;

        approximate_pi_multi(number_of_points);

    return 0;
}
