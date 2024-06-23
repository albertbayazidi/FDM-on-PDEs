#include "visual.h"
#include "utils.h"

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
    // Close file
    fclose(fptr);
}

void saveDDist(int K, float *array, char *filename){

    gsl_matrix *V = gsl_matrix_alloc(2,K);
    gsl_vector_view v = gsl_matrix_row(V,0);

	sigmod_halfs(&v.vector, array); 		

    saveMat(V,filename);

    gsl_matrix_free(V);
}

char* fileloc(char *filename, char *condition){
    // input: filename, folder
    // output: file location REMEMBERED TO FREE MEMORY OF RESULT


    char *result = malloc(100);
    strcpy(result, filename);  // Copy str1 into result
    strcat(result, "_");   // Add a space to result
    strcat(result, condition);  // Concatenate str2 to result
    strcat(result, ".txt");  // Concatenate str2 to result

    return result;
}