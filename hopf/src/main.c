#include "matrix.h"
#include "visual.h"
#include "solve.h"
#include "utils.h"


int main (void){
	int N = 400;  					//number of time steps
	int K = 400;					//number of space steps internal points
	int stop = N-1; 				// number of time steps to compute MAX N-1
	float c = 1;
	
	float gamma = computeGamma(K, N, c);
	float array[4] = {0.5, 10, 1.0, 0.0}; 	// offset, rate, magnitude, minimum

	//solveHopfConst(K, N, gamma, array, stop, roundHat); 
	solveHopfNonConstFlux(K, N, gamma, array, stop, sharpWave);

	return 0;
}

// remember to free vec or mat when not used anymore