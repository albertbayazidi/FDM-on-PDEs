#include "utils.h"

void dirichletDelta(gsl_vector *u,float *array){
    // input: u, x0, sigma
    // u: vector to store the dirichlet delta
    // array: {offset/x0, sigma, K, --} 
    // output: u

    float x0 = array[0];
    float sigma = array[1];
    float K = array[2];

    float fact = 1/(sigma*M_SQRTPI);
    x0 = x0/(K);
    for (float i = 0; i < K; i++)
    {
        float x = i/(K);
        float arg = pow((x-x0)/sigma,2);
        float val = fact*exp(-arg);
        gsl_vector_set(u,i,val);
    }

} 

void initialzeU(gsl_matrix *U,float *array){
    // input: U, N, K, x0, sigma
    // U: matrix to store the solution
    // array = {x0, sigma, --, --}
    // output: U


    gsl_vector_view u0 = gsl_matrix_row(U,0); 

	dirichletDelta(&u0.vector,array);

}

// apply boundary conditions to the matrices (BOTH A AND B)

// neumann
void applayNeumann(double alfa, gsl_matrix *A, gsl_matrix *B){
    // apply neumann boundary conditions to the matrices (BOTH A AND B)
    // input: alfa, A, B

    int K = A->size1;
    gsl_matrix_set(A,0,1,-alfa);
    gsl_matrix_set(A,K-1,K-2,-alfa);

    gsl_matrix_set(B,0,1,alfa);
    gsl_matrix_set(B,K-1,K-2,alfa);
    
}


void applayNeumannNonConstDiff(float beta, gsl_vector *D, gsl_matrix *A, gsl_matrix *B){

    int K = A->size1;
    float D_mZero = gsl_vector_get(D,0)*beta; 
    float D_pZero = gsl_vector_get(D,1)*beta; 
    float D_mLast = gsl_vector_get(D,K-2)*beta;
    float D_pLast = gsl_vector_get(D,K-1)*beta; 

    gsl_matrix_set(A,0,1,-D_mZero - D_pZero);
    gsl_matrix_set(A,K-1,K-2,-D_mLast - D_pLast);
    
    gsl_matrix_set(B,0,1,D_mZero + D_pZero);
    gsl_matrix_set(B,K-1,K-2,D_mLast + D_pLast);

}


// dirichlet
void applayDirichlet(gsl_vector *u){
    // apply dirichlet boundary conditions to the vector u
    // input: u
    
    int n = u->size;
    gsl_vector_set(u,0,0);
    gsl_vector_set(u,n-1,0); // double check this

}

void unboundAnalSol(int N,int K, float D, float *array, char *condition){
    // input: N, D, x0
    // N: number of time steps
    // D: diffusion constant 
    // x0: center of the delta // fixed at centet
    // output: U

    float x0 = array[0];

    gsl_matrix *UAnal = gsl_matrix_alloc(N,K);
    gsl_vector_view u0 = gsl_matrix_row(UAnal,0);
    dirichletDelta(&u0.vector,array);

    x0 = x0/(K);
    for (float i = 1; i < N; i++){
        float fact = 4*D*i/(N);
        for (float j = 0; j < K; j++){
            float x = j/(K);
            float coeff = sqrt(1/(M_PI*fact));
            float arg = pow((x-x0),2)/fact;
            float val = coeff*exp(-arg);
            gsl_matrix_set(UAnal,i,j,val);
        }
    }

    // save the solution
    char filename[] = "misc/data/anal_unbound.txt";

    saveMat(UAnal,filename);

    // free memory
    gsl_matrix_free(UAnal);
}

void sigmod(gsl_vector *D, float *array ){
    // input: D, peak, rate
    // D: vector to store the diffusion 
    // peak: peak of the sigmod
    // rate: the rate at which the function goes from 0 to 1
    // output: deltaD

    float offset = array[0];
    float rate = array[1];
    float magnitude = array[2];
    float minimum = array[3];

    int n = D->size;
    offset = offset/(n);
    for (float i = 0; i < n; i++)
    {
        float x = i/(n);
        float val = magnitude/(1+exp(-rate*(x-offset))) + minimum;
        gsl_vector_set(D,i,val);
    }
}

void sigmod_halfs(gsl_vector *D, float *array ){

    float offset = array[0];
    float rate = array[1];
    float magnitude = array[2];
    float minimum = array[3];

    int n = D->size;
    offset = offset/(n);
    for (float i = 0; i < n; i++)
    {
        float x = (i+0.5)/(n);
        float val = magnitude/(1+exp(-rate*(x-offset))) + minimum;
        gsl_vector_set(D,i,val);
    }

}

void deltaDiff(gsl_matrix *DMat){
    // input: D, N
    // D: vector to store the diffusion 
    // N: number of time steps
    // output: deltaD

    int n = DMat->size2;

    gsl_vector_view DArray = gsl_matrix_row(DMat,0);
	gsl_vector_view deltaDArray = gsl_matrix_row(DMat,1);

    for (int i = 0; i < n-1; i++)
    {
        float val = gsl_vector_get(&DArray.vector,i) + gsl_vector_get(&DArray.vector,i+1); // kan gjøres mer effektiv
        gsl_vector_set(&deltaDArray.vector,i,val);
    }

}

float CFL(int K, int N, float D){
    // input: D, N
    // D: diffusion constant 
    // N: number of time steps
    // output: deltaT

	double deltaT = 1/(double)(N);  //time step
	double deltaX = 1/(double)(K);  //space step

	double alfa = D*deltaT/(pow(deltaX,2));
    return alfa;
}

void pseudoCFL(int K ,int N, float *array, void (*fun)(gsl_vector *, float *)){
    // input: D, N
    // D: vector to store the diffusion 
    // N: number of time steps
    // output: deltaD

    gsl_vector *D = gsl_vector_alloc(K);
    (*fun)(D,array);  
    float Dmax = gsl_vector_max(D);
    
    float deltaT = 1/(float)(N);
    float deltaX = 1/(float)(K);
    
    float val = Dmax*deltaT/pow(deltaX,2);
    printf("Pseudo CFL number: %f \n",val);

    gsl_vector_free(D);

}

float RefBoundConstDX0(int n, int size, float x0){
    // x0 allready scaled
    float vnx0_scale;

    if (n == 0){
        vnx0_scale = 1;              // originaly it shoyld be sqrt(1/L) but L = 1
    }
    else{
        vnx0_scale =  sqrt(2) * cos(n*M_PI*x0); 

    }
    return vnx0_scale;
}

float AbsoBoundConstDX0(int n, int size, float x0){
    // x0 allready scaled
    float vnx0_scale;

    if (n == 0){
        vnx0_scale = 0;              // originaly it shoyld be sqrt(1/L) but L = 1
    }
    else{
        vnx0_scale =  sqrt(2) * sin(n*M_PI*x0); 

    }
    return vnx0_scale;
}

void RefBoundConstD(gsl_vector *x, int n, float x0){
    int size = x->size;
    float vnx0_scale;
    if (n == 0){
        gsl_vector_set_all(x,1);
    }
    else{
        for (float i = 0; i < size; i++){
            gsl_vector_set(x,i,sqrt(2) * cos(n*M_PI*i/(float)size));
        } 
    }

    vnx0_scale = RefBoundConstDX0(n, size, x0);

    gsl_vector_scale(x,vnx0_scale);
}

void AbsoBoundConstD(gsl_vector *x, int n, float x0){
    int size = x->size;
    float vnx0_scale;
    if (n == 0){
        gsl_vector_set_all(x,0);
    }
    else{
        for (float i = 0; i < size; i++){
            gsl_vector_set(x,i,sin(n*M_PI*i/(float)size));
        }

    }
    vnx0_scale =  AbsoBoundConstDX0(n, size, x0);
    gsl_vector_scale(x,vnx0_scale);
}

void boundAnalSolConstD(int N, int K, float D, float *array, char *condition){
    float x0 = array[0];
    x0 = x0/(float)K;

    int n_max = 50;
    gsl_matrix *UAnalBoundConstD = gsl_matrix_alloc(N-1,K); // rem to free after saving
    
    float deltaT = 1/(double)N;
    int i = 0;

    // kanskje ha i egen funk (man kan erstatte RefBoundConstD med en pointer til funksjonen)
    if (strcmp(condition,"neumann") == 0){
        for (float t = deltaT; t <= 1; t = t + deltaT){
            gsl_vector_view sum  = gsl_matrix_row(UAnalBoundConstD,i);
            gsl_vector *v = gsl_vector_alloc(K);  
            for (int n = 0; n < n_max; n++){
                float scale = exp(-pow((n*M_PI),2)*D*t);
                RefBoundConstD(v,n,x0);
                gsl_vector_scale(v,scale);
                gsl_vector_add(&sum.vector,v);
            }
            i++;
            gsl_vector_free(v);
        }
    }
    else{
        for (float t = deltaT; t <= 1; t = t + deltaT){
            gsl_vector_view sum  = gsl_matrix_row(UAnalBoundConstD,i);
            gsl_vector *v = gsl_vector_alloc(K);  
            for (int n = 0; n < n_max; n++){
                float scale = exp(-pow((n*M_PI),2)*D*t);
                AbsoBoundConstD(v,n,x0);
                gsl_vector_scale(v,scale);
                gsl_vector_add(&sum.vector,v);
            }
            i++;
            gsl_vector_free(v);
        }
    }

    // save the solution
    char filename[] = "misc/data/anal_bound_constD";
    char* result = fileloc(filename,condition);
    saveMat(UAnalBoundConstD,result);
    free(result);
    gsl_matrix_free(UAnalBoundConstD);
}


float A_p(float x0,float D_p,float D_m,float t){
    float val = 0;
    float secTerm = gsl_sf_erf(x0/sqrt(4*D_p*t));
    float thirdTerm = sqrt(D_m/D_p)*exp((D_p-D_m)*pow(x0,2)/(4*D_p*D_m*t));
    float trhidScale = (1-gsl_sf_erf(x0/sqrt(4*D_m*t)));

    val = 2/(1+secTerm+thirdTerm*trhidScale);

    //printf("t %f secTerm %f thirdTerm %f trhidScale %f val %f\n",t, secTerm,thirdTerm,trhidScale, val);
    return val;
}

float A_m(float x0,float D_p,float D_m,float t,float A_p){
    float firstScale = sqrt(D_m/D_p)*exp((D_p-D_m)*pow(x0,2)/(4*D_p*D_m*t));
    //printf("t %f firstScale %f A_p %f\n",t, firstScale,A_p);
    float val = A_p*firstScale;
    return val;
}

void unboundAnalSolNonConstD(int N, int K, float *array){
    gsl_vector *D = gsl_vector_alloc(K);
    sigmod(D,array);
    gsl_matrix *U = gsl_matrix_alloc(N,K);
    float x0 = array[0]/K;
    float D_p = gsl_vector_max(D);
    float D_m = gsl_vector_min(D); 

    for (float i = 10; i < N; i++){ //tid 
        gsl_vector_view u = gsl_matrix_row(U,i);
        float t = i/N;
        for (float j = 0; j < K; j++){ //rom
            float x = j/K;
            float a_p = A_p(x0, D_p, D_m, t); 
            
            if (x >= x0){
                float firstScale = sqrt(4*M_PI*D_p*t);
                float secondScale = exp(-pow((x-x0),2)/(4*D_p*t));
                gsl_vector_set(&u.vector,j,a_p/firstScale*secondScale);
            }
            else{
                float a_m = A_m(x0,D_p, D_m, t,a_p);
                float firstScale = sqrt(4*M_PI*D_m*t);
                float secondScale = exp(-pow((x-x0),2)/(4*D_m*t));
                gsl_vector_set(&u.vector,j,a_m/firstScale*secondScale);
            }
        }
    }

    // save the solution
    char filename[] = "misc/data/anal_unbound_nonConstD.txt";
    saveMat(U,filename);

    // freeing memory
    gsl_vector_free(D);
    gsl_matrix_free(U);



}
