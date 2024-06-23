#include "utils.h"


void initialzeU(gsl_matrix *U,float *array,void (*fun)(gsl_vector *,float,float)){
    // input: U, N, K, x0, sigma
    // U: matrix to store the solution
    // array use NULL if you don't need to parametrs {offset,varible}
    // x0: center of the delta
    // sigma: width of the delta should be smaller than 1
    // output: U

	// applay initail conditions on first row of U
	gsl_vector_view u0 = gsl_matrix_row(U,0); 


    int n = u0.vector.size;

    if (array[0] > n || array[0] < 0)
    {
        printf("ERROR: x0 is too big or negative\n");
        exit(1);
    }
    
    (*fun)(&u0.vector,array[0],array[1]);

}

void dirichletDelta(gsl_vector *u,float x0,float sigma){
    // input: u, x0, sigma
    // u: vector to store the dirichlet delta
    // x0: center of the delta
    // sigma: width of the delta should be smaller than 1
    // output: u

    int N = u->size;

    float fact = 1/(sigma*M_SQRTPI);
    x0 = x0/(N);
    for (float i = 0; i < N; i++)
    {
        float x = i/(N);
        float arg = pow((x-x0)/sigma,2);
        float val = fact*exp(-arg);
        gsl_vector_set(u,i,val);
    }

} 

void sinus(gsl_vector *u,float phase,float frequency){
    // input: u, frequency, phase
    // u: vector to store the sinus
    // phase: phase of the sinus (might need to be scaled by N)
    // frequency: frequency of the sinus
    // output: u

    float N = u->size;

    for (int i = 0; i < N; i++){
        float x = i/(N);
        float val = sin(frequency*x + phase);
        gsl_vector_set(u,i,val);
    }

}


void zeroBound(gsl_matrix *U){
    // input: U
    // U: matrix to store the solution
    // output: U
    gsl_matrix_set_zero(U);
}


void initalCondSinSin(gsl_matrix *U,float* array){
    // input: U
    // U: matrix to store the solution
    // output: U

    int N = U->size1 - 1;
    int K = U->size2 - 1;

    for (int i = 1; i < N; i++){ // y = i/
        float y = i/(float)(N);
        for (int j = 1; j < K; j++){  // x = j/
            float x = j/(float)(K);
            gsl_matrix_set(U,i,j,sin(M_PI*x)*sin(2*M_PI*y));
        }
    }
}

void initalCondExp(gsl_matrix *U,float* array){
    // input: U
    // U: matrix to store the solution
    // output: U

    int N = U->size1;
    int K = U->size2;
    float r0_x = array[0];
    float r0_y = array[1];
    float sigma = array[2];

    for (int i = 1; i < N-1; i++){ // y = i/
        float y = i/(float)N;
        for (int j = 1; j < K-1; j++){  // x = j/
            float x = j/(float)K;
            gsl_matrix_set(U,i,j,exp(-(pow(x-r0_x,2)+pow(y-r0_y,2))/sigma));
        }
    }
}

void zeroBoundVec(gsl_vector *v, int K){
    // print a vector as a matrix
    int K2 = v->size; 
    gsl_vector_view firsRow = gsl_vector_subvector(v, 0, K);
    gsl_vector_view lastRow = gsl_vector_subvector(v, K2-K, K);
    gsl_vector_set_zero(&firsRow.vector);
    gsl_vector_set_zero(&lastRow.vector);
    
    unsigned int i = 1;

    while (i<K-1){
        gsl_vector_set(v,i*K,0);
        gsl_vector_set(v,(i)*(K) + (K-1),0);
        i++;
    }
   
}
