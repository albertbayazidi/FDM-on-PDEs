#pragma once
#include <gsl/gsl_linalg.h>

void triDiag(gsl_matrix *A,float scale1, float slcae2, float scale3);

void biDiag(gsl_matrix *B,float subdiagVal,float superdiagVal);
