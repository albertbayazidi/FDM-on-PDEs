#include "matrix.h"
#include "visual.h"
#include "solve.h"
#include "utils.h"


int main (void){

	int N = 10000; 					//number of time steps
	int K = 51;     				//number of space steps internal points
	int stop = N-1;					// number of time steps to compute (max N-1)	
	double D = 1; 				//diffusion constant x²/t 
	char condition[15] = "neumann"; // boundary condition (neumann or dirichlet)

	int x0 = 25;					// center of the delta. location is based vector position
	float sigma = 0.01;			// width of the peak

	float initalCond[4] = {x0,sigma, K,0}; // inital conditions for the delta peak

	// CONSTANT DIFFUSION TASKS #####	
	// gsl_matrix *U = gsl_matrix_alloc(N,K);
	// initialzeU(U,initalCond);
	// constantDiffusionTasks(N, K, D, initalCond, U, condition, stop);
	// boundAnalSolConstD(N, K, D, initalCond, condition);	
	

	// NON CONSTANT DIFFUSION TASKS #####

	int offset = x0;				// center of the peak. Location is based vector position
    float rate = 10;
    float magnitude = 1; 
	float minimum = 0.1;

	float array[4] = {offset,rate,magnitude,minimum};

	pseudoCFL(K, N, array, sigmod);

	char filenamediffDist[] =  "misc/data/diffDist.txt";

	saveDDist(K, array, filenamediffDist);
	gsl_matrix *U = gsl_matrix_alloc(N,K);
	initialzeU(U,initalCond);
	computeSolVariableDiff(N, K, array, U, condition, stop);	

	unboundAnalSolNonConstD(N, K, array);

	return 0;
}
