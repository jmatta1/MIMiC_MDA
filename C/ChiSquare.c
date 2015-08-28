//uncomment these for lots of error messages and messages about exact status
#define MAKE_AND_FREE_OUTPUT


#include"./ChiSquare.h"

#include<stdio.h>
#include<stdlib.h>


void* makeMdaStruct(int numPts, int numL)
{
    MdaData* temp = (MdaData*) malloc(sizeof(MdaData));
#ifdef MAKE_AND_FREE_OUTPUT
    printf("MdaDataStruct is: %d bytes\n", sizeof(MdaData));
    printf("MdaDataStruct is at: %p\n", temp);
#endif
    temp->numPoints = numPts;
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Number of points is: %d\n", numPts);
#endif
    temp->numLVals = numL;
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Number of L values is: %d\n", numL);
#endif
    temp->data = (float*) malloc(sizeof(float)*numPts);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Data size is: %d byes\n", sizeof(float)*numPts);
    printf("Data is at: %p\n", temp->data);
#endif
    temp->distArrays = (float*) malloc(sizeof(float)*numPts*numL);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Distribution size is: %d byes\n", sizeof(float)*numPts*numL);
    printf("Distributions are at: %p\n", temp->distArrays);
#endif
    return ((void*)temp);
}

void freeMdaStruct(void* strPtr)
{
    MdaData* temp = (MdaData*)strPtr;
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the distributions\n");
#endif
    free(temp->distArrays);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the data\n");
#endif
    free(temp->data);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the structure\n");
#endif
    free(temp);
}

