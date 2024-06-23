#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "visual.h"


void dirichletDelta(gsl_vector *u,float *array); // when you have time and are intrested optimize this to take parameters as array sa CFL

void sigmod(gsl_vector *D,float *array);

void sigmod_halfs(gsl_vector *D, float *array );

void deltaDiff(gsl_matrix *DMat);

void initialzeU(gsl_matrix *U,float *array);

void applayNeumann(double alfa,gsl_matrix *A, gsl_matrix *B);

void applayNeumannNonConstDiff(float beta, gsl_vector *D, gsl_matrix *A, gsl_matrix *B);

void applayDirichlet(gsl_vector *u);

void unboundAnalSol(int N,int K, float D, float* array, char *condition);

float CFL(int K, int N, float D);

void pseudoCFL(int K, int N, float *array, void (*fun)(gsl_vector *, float *));

float AbsoBoundConstDX0(int n, int size, float x0);

float RefBoundConstDX0(int n, int size, float x0);

void RefBoundConstD(gsl_vector *x, int n, float x0);

void AbsoBoundConstD(gsl_vector *x, int n, float x0);

void boundAnalSolConstD(int N, int K, float D, float *array, char *condition);

void unboundAnalSolNonConstD(int N, int K, float *array);

float A_p(float x0,float D_p,float D_m,float t);

float A_m(float x0,float D_p,float D_m,float t,float A_p);