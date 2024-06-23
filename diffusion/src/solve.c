#include "solve.h"
#include "visual.h"
#include "utils.h"

void constantDiffusionTasks(int N, int K, float D, float *array ,gsl_matrix *U, char *condition, int stop){
    // solve equation with constant diffusion constant

    // compute the solution for constant diffusion
    float alfa = CFL(K, N, D);
    printf("CFL number: %f \n",alfa);
	computeSolConstDiff(alfa, N, K, U, condition, stop);

    // save the solution
    char filenameCont[] = "misc/data/numsol";
    char* result = fileloc(filenameCont, condition);
	saveMat(U,result);

    // free memory
    free(result);
    gsl_matrix_free(U);

    // compute the analytical solution for unbound case used for task 2.7
	unboundAnalSol(N,K,D,array,condition);

    // should make analytical solution for bound case both dirichlet and neumann
    /////////////

}

void computeSolConstDiff(double alfa, int N, int K, gsl_matrix *U, char *condition ,int stop){
    // solve equation with constant diffusion constant
    // input: alfa, N, K, U
    // alfa: scaled diffusion constant 
    // N: number of time steps
    // K: number of space steps
    // U: matrix to store solution with inital conditions applied on first rowz
    // output: U
    
	gsl_matrix *A = gsl_matrix_alloc(K,K);          // LHS matrix
    gsl_matrix *B = gsl_matrix_alloc(K,K);          // RHS matrix

    initialzeConstA(alfa, A);                       // initialize the diagonals of A
    initialzeConstB(alfa, B);                       // initialize the diagonals of B
    

    // applay the boundary conditions (on both A AND B)
    if (strcmp(condition,"neumann") == 0)
    {
        applayNeumann(alfa,A,B);
    }

    gsl_vector_view diagA = gsl_matrix_diagonal(A);
    gsl_vector_view subDiagA = gsl_matrix_subdiagonal(A,1);
    gsl_vector_view superDiagA = gsl_matrix_superdiagonal(A,1);

    for (int j = 0; j < stop; j++){ 
    	gsl_vector_view u = gsl_matrix_row(U,j); 
		gsl_vector_view uNext = gsl_matrix_row(U,j+1);

        if (strcmp(condition,"dirichlet") == 0){
            applayDirichlet(&u.vector); // under construction
        }
        // compute RHS called v
		makeCrankNicolsonRHS(B,&u.vector,&uNext.vector);
		
		// solve system of equations Au = v
		gsl_linalg_solve_tridiag(&diagA.vector,&superDiagA.vector,&subDiagA.vector,&uNext.vector,&uNext.vector);

    }

    gsl_matrix_free(A);
    gsl_matrix_free(B);
}




void computeSolVariableDiff(int N, int K,float *array, gsl_matrix *U, char *condition, int stop){
    // solve equation with variable diffusion constant
    // input: D, N, K, U
    char filenameNonConst[] = "misc/data/numsolNonConst";

    double deltaT = 1/(double)(N);              //time step
	double deltaX = 1/(double)(K);              //space step

	double beta = deltaT/(2*pow(deltaX,2));


    gsl_matrix *A = gsl_matrix_alloc(K,K);
    gsl_matrix *B = gsl_matrix_alloc(K,K);


	gsl_matrix *DMat = gsl_matrix_alloc(2,K+1); // first row is D_{i+0.5} second row is deltaD
	gsl_vector_view DArray = gsl_matrix_row(DMat,0);
	
	sigmod_halfs(&DArray.vector, array); 		// pick distribution of diffusion

	deltaDiff(DMat);

	diffMat(DMat, beta, A, "A");
	diffMat(DMat, beta, B, "B");

	gsl_vector_view diagA = gsl_matrix_diagonal(A);
    gsl_vector_view subDiagA = gsl_matrix_subdiagonal(A,1);
    gsl_vector_view superDiagA = gsl_matrix_superdiagonal(A,1);

	applayNeumannNonConstDiff(beta, &DArray.vector, A, B);
	for (int j = 0; j < stop; j++){ 
    	gsl_vector_view u = gsl_matrix_row(U,j); 
		gsl_vector_view uNext = gsl_matrix_row(U,j+1);

        if (strcmp(condition,"dirichlet") == 0){
            applayDirichlet(&u.vector); // under construction
        }
        // compute RHS called v
		makeCrankNicolsonRHS(B,&u.vector,&uNext.vector);
		
		// solve system of equations Au = v
		gsl_linalg_solve_tridiag(&diagA.vector,&superDiagA.vector,&subDiagA.vector,&uNext.vector,&uNext.vector);
    }
    
    // saving result
    char* result = fileloc(filenameNonConst, condition);
	saveMat(U,result);

    // freeing memory
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    free(result);

}