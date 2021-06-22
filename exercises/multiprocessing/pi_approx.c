#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "random_number.h"
#include "pi_approx.h"

#define NUMBER_OF_THREADS 2

typedef struct place_points_struct
	{
	unsigned int seed;
	double number_of_points;
	int *points_in_circle;
} place_points_struct;

void *place_points(void *place_points_struct_input){
	place_points_struct *arguments = (place_points_struct *) place_points_struct_input;
	unsigned int seed = arguments->seed;
	double number_of_points = arguments->number_of_points;
	int *points_in_circle = arguments->points_in_circle;
	double circle_radius = 1.0;
	double x=0;
	double y=0;

for (int i = 0; i < number_of_points; i++)
    {
        //Generate random points
        x = random_number(&seed);
        y = random_number(&seed);

        if (sqrt(pow(x, 2) + pow(y, 2)) <= circle_radius)
        {
            (*points_in_circle)++;
        }
    }
    pthread_exit((void *) place_points_struct_input);
    return NULL;
}



void approximate_pi_multi(double number_of_points)
{
    pthread_t thread_array[NUMBER_OF_THREADS];
    pthread_attr_t attributes;
    place_points_struct *place_points_struct_array = malloc(NUMBER_OF_THREADS * sizeof(place_points_struct));
    int error;
    void *status;

    //Initialize threads
    pthread_attr_init(&attributes);
    pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);

    unsigned int master_seed = time(NULL);
    int *points_in_circle_array = malloc(NUMBER_OF_THREADS * sizeof(int));

    //Create threads that run in parallel in loop
    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        unsigned int thread_seed = master_seed + i; //Set unique seed for each thread
        points_in_circle_array[i] = 0;

        place_points_struct_array[i].seed = thread_seed;
        place_points_struct_array[i].number_of_points = number_of_points;
        place_points_struct_array[i].points_in_circle = &points_in_circle_array[i];

        error = pthread_create(&thread_array[i], &attributes, place_points, (void *) &(place_points_struct_array[i]));
        if (error)
        {
            printf("pthread_join() returned failure");
            exit(-1);
        }

    }
    pthread_attr_destroy(&attributes);

    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        error = pthread_join(thread_array[i], &status);
        if (error)
        {
            printf("pthread_join() returned failure");
            exit(-1);
        }

	}
}
