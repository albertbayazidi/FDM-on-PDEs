#pragma once
#include <gsl/gsl_linalg.h>
#include <string.h>
#include "visual.h"


void makeTriMat(gsl_vector *diag,gsl_vector *subDiag,gsl_vector *superDiag,gsl_matrix *A);

void makeCrankNicolsonRHS(gsl_matrix *B, gsl_vector *u, gsl_vector *v);

void initialzeConstA(double alfa, gsl_matrix *A);

void initialzeConstB(double alfa, gsl_matrix *B);


void diffMat(gsl_matrix *DMat, float beta, gsl_matrix *Mat,char *condition);



