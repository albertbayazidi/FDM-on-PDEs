#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

void roundHat(gsl_vector *v, float* array);

void sigmod(gsl_vector *D, float *array);

void waveFront(gsl_vector *v, float* array);

void sharpWave(gsl_vector *v, float* array);
