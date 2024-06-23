#include <stdio.h>
#include "visual.h"


void printVec(gsl_vector *v){
    // print a vector
    for (unsigned int i = 0; i < v->size; i++)
    {
        printf("[%d], %f \n",i,gsl_vector_get(v,i));
    }
    
}

void printMat(gsl_matrix *A){
    // print a matrix
    for (unsigned int i = 0; i < A->size1; i++)
    {   
        printf("[%d], ",i);
        for (unsigned int j = 0; j < A->size2; j++)
        {
            printf("%f ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
}


void saveMat(gsl_matrix *A, char *filename){
    // save a matrix to a file as a array (driffrent from to other folders)
    // MÅ påasse på at alle elementete endrer seg da A er samme matrise hele tiden
    // det betyr at vi hvis 1,1 ikke endrer seg fra u0 til t1 blir den husket (går bra i teorien)
    
    FILE *fptr;
    int size = A->size1*A->size2;
    // Open file
    fptr = fopen(filename, "a");

    // Write 
    for (unsigned int i = 0; i < size; i++){   
        fprintf(fptr,"%f,",A->data[i]);
    }
    fprintf(fptr, "\n");

    // Close file
    fclose(fptr);
}

void printVecAsMat(gsl_vector *v, int K){
    // print a vector as a matrix
    for (unsigned int i = 0; i < v->size; i++)
    {
        printf("%f ",gsl_vector_get(v,i));
        if ((i+1)%K == 0)
        {
            printf("\n");
        }
    }
    printf("\n");
}
