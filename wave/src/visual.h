#pragma once
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <stdlib.h>

#define SHELLSCRIPT "\
#/bin/bash \n\
rm \"misc/data/numSolSinInit.txt\" \n\
"


void printVec(gsl_vector *v);

void printMat(gsl_matrix *A);

void saveMat(gsl_matrix *A, char *filename);

void printVecAsMat(gsl_vector *v, int K);