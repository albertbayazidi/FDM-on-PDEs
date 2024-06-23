#include "matrix.h"
#include "visual.h"
#include "solve.h"
#include "utils.h"

int main (void){
	int K = 20;  				// number of x steps
	int T = 50;					// number of time steps
	float c = 1; 				// wave speed
	float array[3] = {0.5,0.5,0.001};

    system(SHELLSCRIPT);		// remove privouse numsol

	solveWaveEq(K, T, c, array, initalCondSinSin);

	return 0;
}


