#include "matrix.h"

void triDiag(gsl_matrix *A,float diagVal,float subdiagVal, float superdiagVal){
    // input: A, scales 
    // A: matrix to store the tridiagonal matrix
    // scales: the scales of the terms in the equation
    gsl_vector_view diag = gsl_matrix_diagonal(A);
    gsl_vector_view subdiag = gsl_matrix_subdiagonal(A,1);
    gsl_vector_view superdiag = gsl_matrix_superdiagonal(A,1);

    gsl_vector_set_all(&diag.vector,diagVal);
    gsl_vector_set_all(&subdiag.vector,subdiagVal);
    gsl_vector_set_all(&superdiag.vector,superdiagVal);
}

void biDiag(gsl_matrix *B,float subdiagVal,float superdiagVal){
    // input: A, scales 
    // A: matrix to store the tridiagonal matrix
    // scales: the scales of the terms in the equation
    gsl_vector_view subdiag = gsl_matrix_subdiagonal(B,1);
    gsl_vector_view superdiag = gsl_matrix_superdiagonal(B,1);

    gsl_vector_set_all(&subdiag.vector,subdiagVal);
    gsl_vector_set_all(&superdiag.vector,superdiagVal);
}