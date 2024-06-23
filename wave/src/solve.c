#include "solve.h"


float computeBeta(int K, int T, float c){    
    float h = 1.0/(float)(K-1); // space step
    float dt = 1.0/(float)(T-1); // time step
    float beta = pow(c,2)*pow(dt,2)/(float)pow(h,2); // courant number

    printf("beta = %f,h = %f, dt = %f \n",beta,h,dt);
	return beta;
}

void solveWaveEq(int K, int T, float c,float *array,void (*fun)(gsl_matrix*, float*)){
	float N = K;											// number of space_steps x	
    float K2 = pow(K,2);									// number of space_steps y 	
	float beta = computeBeta(K, T, c);						// compute beta
    
    /// inisialize U data
	gsl_matrix *U = gsl_matrix_alloc(N,K); 					// matrix to store the solution 
	zeroBound(U);
	(*fun)(U,array);
	saveMat(U,"misc/data/numSolSinInit.txt");	

	gsl_matrix *A = gsl_matrix_alloc(K2,K2); 			

	// inisialize matrix
	pentaDiag(A,beta,K); 

	// convert U to vector u0
	gsl_vector_view u0 = gsl_vector_view_array(U->data, U->size1 * U->size2); 

	// compute uMinus
	gsl_vector *uMinus = gsl_vector_alloc(K2); 				
	gsl_blas_daxpy(-1,&u0.vector,uMinus); 					// uMinus = -u0 //
	gsl_blas_dgemv(CblasNoTrans,1,A,&u0.vector,1.0,uMinus); //
	zeroBoundVec(uMinus, K);								// setter boundry 0

	// solve the wave equation and frees memory
	solveLoop(T, U, A, uMinus);	
}

void solveLoop(int T, gsl_matrix *U, gsl_matrix *A, gsl_vector *uMinus){
	// solve the wave equation
	// input: K, U, A, I, uMinus, uPrev, uCur, uNew, u0
	// also saves the solution to file and frees memory
	int count = 0;
	int K = U->size1; 										// number of space_steps x 
	int K2 = uMinus->size;

	// convert matrix U to vectors
	gsl_vector *uPrev = gsl_vector_alloc(K2);
	gsl_vector *uCur = gsl_vector_alloc(K2); 
	gsl_vector_memcpy(uPrev,uMinus);									// uPrev = uNew	

	while (count < T-1){ 
		gsl_vector_view uNew = gsl_vector_view_array(U->data, U->size1 * U->size2);
		gsl_vector_memcpy(uCur,&uNew.vector);
		gsl_vector_set_zero(&uNew.vector);								// uNew = 0

		gsl_blas_daxpy(-1,uPrev,&uNew.vector);							// uNew = -uPrev KANSKJE DETTE MÅ LAGRES I NOE ANNET ENN EN VEK VIEW
		gsl_blas_dgemv(CblasNoTrans,1,A,uCur,1.0,&uNew.vector);			// uNew = uNew + A*uCur

		zeroBoundVec(&uNew.vector, K);
		gsl_vector_memcpy(uPrev, uCur);							// uPrev = uCur
		gsl_vector_memcpy(uCur,&uNew.vector);
		
		saveMat(U,"misc/data/numSolSinInit.txt");
		count++;
	}

	gsl_vector_free(uCur); 				// vectors
	gsl_vector_free(uPrev);				// vectors
	gsl_vector_free(uMinus);
	
	gsl_matrix_free(U);					// matrices
	gsl_matrix_free(A);					// matrices
}