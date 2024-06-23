#pragma once
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "matrix.h"
#include "visual.h"
#include "utils.h"

float computeBeta(int K, int T, float c);


void solveWaveEq(int K, int T, float c, float *array,void (*fun)(gsl_matrix*, float*));


void solveLoop(int T, gsl_matrix *U, gsl_matrix *A, gsl_vector *uMinus);