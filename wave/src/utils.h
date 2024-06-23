#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include "visual.h"
#include "matrix.h"

void dirichletDelta(gsl_vector *u,float x0,float sigma);

void sinus(gsl_vector *u,float frequency,float phase);

void initialzeU(gsl_matrix *U,float *array,void (*fun)(gsl_vector *,float ,float));

void zeroBound(gsl_matrix *U);

void zeroBoundVec(gsl_vector *u, int K);

void initalCondSinSin(gsl_matrix *U, float* array);

void initalCondExp(gsl_matrix *U, float* array);