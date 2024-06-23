#include "matrix.h"

void pentaDiag(gsl_matrix *A,float beta, int K){
    // makes simple tridiagonal matrix with -2 on diagonal and 1 on sub and super diagonal

    gsl_vector_view diag = gsl_matrix_diagonal(A);
    gsl_vector_view subdiag = gsl_matrix_subdiagonal(A,1);
    gsl_vector_view superdiag = gsl_matrix_superdiagonal(A,1);
    gsl_vector_view subNdiag = gsl_matrix_subdiagonal(A,K);
    gsl_vector_view superNdiag = gsl_matrix_superdiagonal(A,K);

    gsl_vector_set_all(&diag.vector,2*(1-2*beta));
    gsl_vector_set_all(&subdiag.vector,beta);
    gsl_vector_set_all(&superdiag.vector,beta);
    gsl_vector_set_all(&subNdiag.vector,beta);
    gsl_vector_set_all(&superNdiag.vector,beta);

}

