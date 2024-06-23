#pragma once
#include <stdio.h>
#include <gsl/gsl_linalg.h>

void printVec(gsl_vector *v);

void printMat(gsl_matrix *A);

void saveMat(gsl_matrix *A, char *filename);