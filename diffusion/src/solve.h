#pragma once
#include <gsl/gsl_linalg.h>
#include <string.h>
#include "matrix.h"

void computeSolConstDiff(double alfa, int N, int K, gsl_matrix *U,char *condition ,int stop);

void constantDiffusionTasks(int N, int K, float D, float *array, gsl_matrix *U, char *condition, int stop);

void computeSolVariableDiff( int N, int K, float *array, gsl_matrix *U, char *condition, int stop);

