#include "matrix.h"

void makeCrankNicolsonRHS(gsl_matrix *B, gsl_vector *u, gsl_vector *v){
    // make the right hand side of the Crank Nicolson equation Bu(n) = v(n)

    gsl_blas_dgemv(CblasNoTrans,1,B,u,0,v);

}

void initialzeConstA(double alfa, gsl_matrix *A){
    // initialize the vectors diagA and subDiagA 

    gsl_vector_view diagA = gsl_matrix_diagonal(A);
	gsl_vector_view subDiagA = gsl_matrix_subdiagonal(A,1);
	gsl_vector_view superDiagA = gsl_matrix_superdiagonal(A,1);

	gsl_vector_set_all(&diagA.vector,alfa+1);
	gsl_vector_set_all(&subDiagA.vector,-alfa/2);
	gsl_vector_set_all(&superDiagA.vector,-alfa/2);

}

void initialzeConstB(double alfa, gsl_matrix *B){
    // initialize the matrices B

    gsl_vector_view diagB = gsl_matrix_diagonal(B);
	gsl_vector_view subDiagB = gsl_matrix_subdiagonal(B,1);
	gsl_vector_view superDiagB = gsl_matrix_superdiagonal(B,1);

    gsl_vector_set_all(&diagB.vector,1-alfa);
	gsl_vector_set_all(&subDiagB.vector,alfa/2);
	gsl_vector_set_all(&superDiagB.vector,alfa/2);

}


void diffMat(gsl_matrix *DMat, float beta, gsl_matrix *Mat, char *condition){
	// initialize the matrix A for variable diffusion
	// DMat: matrix with the diagonal elements of A
	// A: matrix to initialize
	// beta: scaling factor for the diffusion constant beta = deltaT/2(deltaX^2)

	int K = Mat->size1;

	gsl_vector_view DArray = gsl_matrix_row(DMat,0);
	gsl_vector_view deltaDArray = gsl_matrix_subrow(DMat,1,0,K);
	
	if (strcmp(condition,"B") == 0){
		beta = -beta;
	}

	gsl_vector_view diagMat = gsl_matrix_diagonal(Mat);
	gsl_vector_view subDiagMat = gsl_matrix_subdiagonal(Mat,1);
	gsl_vector_view superDiagMat = gsl_matrix_superdiagonal(Mat,1);

	gsl_vector_view subSuperDiagVals = gsl_vector_subvector(&DArray.vector,1,K-1);

	
	gsl_vector_memcpy(&subDiagMat.vector,&subSuperDiagVals.vector); 
	gsl_vector_memcpy(&superDiagMat.vector,&subSuperDiagVals.vector);
	gsl_vector_memcpy(&diagMat.vector,&deltaDArray.vector);	


	gsl_vector_scale(&subDiagMat.vector,-beta);
	gsl_vector_scale(&superDiagMat.vector,-beta);
	gsl_vector_scale(&diagMat.vector,beta);
	gsl_vector_add_constant(&diagMat.vector,1);


}

