#include "solve.h"
#include "visual.h"

float computeGamma(int K, int N, float c){
	// Compute the Courant number
	// input: K, N, c
	// K: number of space steps internal points
	// N: number of time steps
	// c: speed
	// output: gamma

	float deltaT = 1/(float)(N); //time step
	float deltaX = 1/(float)(K); //space step
	float gamma = c*deltaT/deltaX;
	printf("gamma: %f\n",gamma);
	return gamma;
}

void solveHopfConst(int K, int N, float gamma, float* array, int stop, void (*fun)(gsl_vector *, float*)){
    // Solver
	// input: U, K, scaleFirst, scaleSecond, stop
	// U: matrix to store the solution
	// K: number of space steps internal points 
	// stop: number of time steps to compute
	// output: U

	gsl_matrix *U = gsl_matrix_alloc(N,K);
	gsl_matrix *A = gsl_matrix_alloc(K,K);

	// initial condition
	gsl_vector_view u0 = gsl_matrix_row(U,0);
	fun(&u0.vector,array);


	float diagVal = -pow(gamma,2)+1;
	float subdiagVal = gamma/(float)2*(gamma+1);
	float superdiagVal = gamma/(float)2*(gamma-1);

	triDiag(A,diagVal,subdiagVal,superdiagVal);

	for (int j=0; j < stop; j++){
		gsl_vector_view u = gsl_matrix_row(U,j); 
		gsl_vector_view uNext = gsl_matrix_row(U,j+1);

		gsl_blas_dgemv(CblasNoTrans,1.0,A,&u.vector,0.0,&uNext.vector);
		gsl_vector_set(&uNext.vector,0,gsl_vector_get(&u.vector,K-1));
	}

	saveMat(U,"misc/data/numsolConstFlux.txt");
	gsl_matrix_free(A);
	gsl_matrix_free(U);    
}


void solveHopfNonConstFlux(int K, int N, float gamma, float* array, int stop, void (*fun)(gsl_vector *, float*)){
	// rem to free 
	// c must be 1 
	gsl_matrix *U = gsl_matrix_alloc(N,K);
	gsl_matrix *A = gsl_matrix_alloc(K,K);
	gsl_matrix *B = gsl_matrix_alloc(K,K);
	gsl_matrix *BHat = gsl_matrix_alloc(K,K);

	float sigma =  gamma/(float)4;
	printf("sigma: %f\n",sigma);
	// initial condition
	gsl_vector_view u0 = gsl_matrix_row(U,0);
	fun(&u0.vector,array);

	// initialaize A and B
	triDiag(A,-2,1,1);
	biDiag(B,1,1);
	biDiag(BHat,-1,1);
	
	// memory is freed in solveLoopNonConstFlux
	solveLoopNonConstFlux(K, stop, U, A, B, BHat, sigma);
}

void solveLoopNonConstFlux(int K, int stop, gsl_matrix *U, gsl_matrix *A, gsl_matrix *B, gsl_matrix *BHat, float sigma){
	// dummy vector

	gsl_vector *upow = gsl_vector_alloc(K);
	gsl_vector *vec = gsl_vector_alloc(K);
	for (unsigned int j = 0; j < stop; j ++){

		gsl_vector_view u = gsl_matrix_row(U,j); 
		gsl_vector_view uNext = gsl_matrix_row(U,j+1);

		gsl_vector_memcpy(upow,&u.vector);					// upow = u 
		gsl_vector_mul(upow,&u.vector);  					//upow = u²

		gsl_blas_daxpy(1.0,&u.vector,&uNext.vector);		// uNext = u
		
		gsl_blas_dgemv(CblasNoTrans,-sigma,BHat,upow,1,&uNext.vector); // uNext = uNext - sigma*BHat*u²

		gsl_blas_dgemv(CblasNoTrans,-2*pow(sigma,2),B,&u.vector,0.0,vec); // vec = -2*sigma²*B*u
		gsl_vector_mul(vec,upow); 							// vec = -2*sigma²*B*u⊗u² 
		gsl_blas_daxpy(1.0,&uNext.vector,vec);				// uNext = uNext - 2*sigma²*B*u⊗u²

		gsl_vector_set_zero(vec);
		gsl_blas_dgemv(CblasNoTrans,2*pow(sigma,2),B,upow,0.0,vec); // vec = 2*sigma²*B*u²
		gsl_vector_mul(vec,&u.vector);						// vec = 2*sigma²*B*u²⊗u
		gsl_blas_daxpy(1.0,&uNext.vector,vec);				// uNext = uNext + 2*sigma²*B*u²⊗u

		gsl_vector_mul(upow,&u.vector);						// upow = u³
		gsl_blas_dgemv(CblasNoTrans,2*pow(sigma,2),A,upow,1.0,&uNext.vector); // uNext = uNext - 2*sigma²*A*u³
	}

	saveMat(U,"misc/data/numsolNonConstFluxSharpWave.txt");

	// free memory
	gsl_matrix_free(A);
	gsl_matrix_free(U);    
	gsl_matrix_free(B);
	gsl_matrix_free(BHat);
	gsl_vector_free(upow);
	gsl_vector_free(vec);
}