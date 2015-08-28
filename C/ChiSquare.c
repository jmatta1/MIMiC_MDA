// uncomment these for lots of diagnostic messages
//#define MAKE_AND_FREE_OUTPUT
//#define ASSIGN_DATA_OUTPUT

#include"./ChiSquare.h"
#if defined(MAKE_AND_FREE_OUTPUT) || defined(ASSIGN_DATA_OUTPUT)
#include<stdio.h>
#endif
#include<stdlib.h>


void* makeMdaStruct(int numPts, int numL)
{
    MdaData* temp = (MdaData*) malloc(sizeof(MdaData));
    temp->numPts = numPts;
    temp->numLs = numL;
    temp->data = (float*) malloc(sizeof(float)*numPts);
    temp->dists = (float*) malloc(sizeof(float)*numPts*numL);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("MdaDataStruct is: %d bytes\n", sizeof(MdaData));
    printf("MdaDataStruct is at: %p\n", temp);
    printf("Number of points is: %d\n", numPts);
    printf("Number of L values is: %d\n", numL);
    printf("Data size is: %d byes\n", sizeof(float)*numPts);
    printf("Data is at: %p\n", temp->data);
    printf("Distribution size is: %d byes\n", sizeof(float)*numPts*numL);
    printf("Distributions are at: %p\n", temp->dists);
#endif
    return ((void*)temp);
}

void freeMdaStruct(void* strPtr)
{
    MdaData* temp = (MdaData*)strPtr;
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the distributions\n");
#endif
    free(temp->dists);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the data\n");
#endif
    free(temp->data);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the structure\n");
#endif
    free(temp);
}

void setData(void* strPtr, float* divData)
{
    MdaData* data = (MdaData*)strPtr;
#ifdef ASSIGN_DATA_OUTPUT
    printf("struct ptr is: %p\n", data);
    printf("data ptr is: %p\n", data->data);
#endif
    for (int i=0; i<data->numPts; ++i)
    {
        data->data[i] = divData[i];
#ifdef ASSIGN_DATA_OUTPUT
        printf("At cell %d value: %f\n", i, data->data[i]);
        printf("Orig value: %f\n", divData[i]);
#endif
    }
}


