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

    FILE *fptr;

    // Open file
    fptr = fopen(filename, "w");

    // Write 

        for (unsigned int i = 0; i < A->size1; i++)
    {   
        for (unsigned int j = 0; j < A->size2; j++)
        {
            fprintf(fptr,"%f,", gsl_matrix_get(A,i,j));
        }
        fprintf(fptr, "\n");

    }
    // close file
    fclose(fptr);
}

