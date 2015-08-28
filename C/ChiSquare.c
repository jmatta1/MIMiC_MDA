//uncomment these for lots of error messages and messages about exact status
#define INIT_AND_FREE_OUTPUT


#include"./ChiSquare.h"

#include<stdio.h>
#include<stdlib.h>


void* initMdaStruct(int numPts, int numL)
{
    MdaDataStruct* temp = (MdaDataStruct*) malloc(sizeof(MdaDataStruct));
#ifdef INIT_AND_FREE_OUTPUT
    printf("MdaDataStruct is: %d bytes", sizeof(MdaDataStruct));
    printf("MdaDataStruct is at: %d", temp);
#endif
    temp->numPoints = numPts;
#ifdef INIT_AND_FREE_OUTPUT
    printf("Number of points is: %d", numPts);
#endif
    temp->numLVals = numL;
#ifdef INIT_AND_FREE_OUTPUT
    printf("Number of L values is: %d", numL);
#endif
    temp->data = (float*) malloc(sizeof(float)*numPts);
#ifdef INIT_AND_FREE_OUTPUT
    printf("Data size is: %d byes", sizeof(float)*numPts);
    printf("Data is at: %d", (int)temp->data);
#endif
    temp->distArrays = (float*) malloc(sizeof(float)*numPts*numL);
#ifdef INIT_AND_FREE_OUTPUT
    printf("Distribution size is: %d byes", sizeof(float)*numPts*numL);
    printf("Distributions are at: %d", (int)temp->distArrays);
#endif
    return ((void*)temp);
}

void freeMdaStruct(void* strPtr)
{
    MdaDataStruct* temp = (MdaDataStruct*)strPtr;
#ifdef INIT_AND_FREE_OUTPUT
    printf("Freeing the distributions");
#endif
    free(temp->distArrays);
#ifdef INIT_AND_FREE_OUTPUT
    printf("Freeing the data");
#endif
    free(temp->data);
#ifdef INIT_AND_FREE_OUTPUT
    printf("Freeing the structure");
#endif
    free(temp);
}

