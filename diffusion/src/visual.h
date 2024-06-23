#pragma once
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_linalg.h>

char* fileloc(char *filename, char *condition);

void printVec(gsl_vector *v);

void printMat(gsl_matrix *A);

void saveMat(gsl_matrix *A, char *filename);

void saveDDist(int K, float *array, char *filename);