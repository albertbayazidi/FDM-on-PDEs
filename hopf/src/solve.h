#pragma once
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "matrix.h"

float computeGamma(int K, int N, float c);

void solveHopfConst(int K, int N,float gamma, float* array, 
                    int stop, void (*fun)(gsl_vector *, float* ));


void solveHopfNonConstFlux(int K, int N, float gamma, float* array, int stop, void (*fun)(gsl_vector *, float*));

void solveLoopNonConstFlux(int K, int stop, gsl_matrix *U, gsl_matrix *A, gsl_matrix *B, gsl_matrix *BHat, float sigma);
