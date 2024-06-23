#include "utils.h"
#include "visual.h"

void roundHat(gsl_vector *v, float* array){
    // input: v, array
    // v: vector to store the initial condition
    // array: array with the initial condition

    int K = v->size - 1;
    float offset = array[0];

    for (unsigned int i = 0; i <= K; i++){
        float x = i/(float)K;
        float val = 1-pow(4*(x-offset),2);
        if (val < 0){
            val = 0;
        }
        gsl_vector_set(v,i,val);
    }
}
void waveFront(gsl_vector *v, float* array){
    // input: v, array
    // v: vector to store the initial condition
    // array: array with the initial condition

    int K = v->size - 1;
    for (unsigned int i = 0; i <= K; i++){
        float x = i/(float)K;
        if (x < 0.5){
            float val = 1/(1+exp(-30*(x-0.2)));
            gsl_vector_set(v,i,val);
        }
        else{
            float val = 1/(1+exp(30*(x-0.8)));
            gsl_vector_set(v,i,val);
        }
    }
}

void sharpWave(gsl_vector *v, float* array){
    int K = v->size -1;
    for (unsigned int i = 0; i <= K; i++){
        float x = i/(float)K;
        float val = 1/(exp(100*pow((x-0.5),2)));
        gsl_vector_set(v,i,val);
        
    }
}

void sigmod(gsl_vector *v, float *array){
    // input: D, peak, rate
    // D: vector to store the diffusion 
    // peak: peak of the sigmod
    // rate: the rate at which the function goes from 0 to 1
    // output: deltaD

    float offset = array[0];
    float rate = array[1];
    float magnitude = array[2];
    float minimum = array[3];

    int n = v->size-1;
    for (float i = 0; i <= n; i++)
    {
        float x = i/(n);
        float val = magnitude/(1+exp(-rate*(x-offset))) + minimum;
        gsl_vector_set(v,i,val);
    }
}