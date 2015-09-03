// uncomment these for lots of diagnostic messages
//#define MAKE_AND_FREE_OUTPUT
//#define ASSIGN_DATA_OUTPUT
//#define ASSIGN_DIST_OUTPUT
//#define CALC_CHI_OUTPUT

#include"./ChiSquare.h"
#if defined(MAKE_AND_FREE_OUTPUT) || defined(ASSIGN_DATA_OUTPUT) || \
    defined(ASSIGN_DIST_OUTPUT)
#include<stdio.h>
#endif
//removing stdlib since instead of malloc we are now using memalign
//#include<stdlib.h>
#include<malloc.h>

#define ALIGN_BOUND 64

__attribute__((malloc)) void* makeMdaStruct(int numPts, int numL)
{
    MdaData* temp = (MdaData*) malloc(sizeof(MdaData));
    temp->numPts = numPts;
    temp->numLs = numL;
    //temp->data = (float*) malloc(sizeof(float)*numPts);
    //temp->dists = (float*) malloc(sizeof(float)*numPts*numL);
    //temp->resids = (float*) malloc(sizeof(float)*numPts);
    temp->data = (float*) memalign(ALIGN_BOUND,sizeof(float)*numPts);
    temp->dists = (float*) memalign(ALIGN_BOUND,sizeof(float)*numPts*numL);
    temp->resids = (float*) memalign(ALIGN_BOUND,sizeof(float)*numPts);
#ifdef MAKE_AND_FREE_OUTPUT
    printf("\n\nMdaDataStruct is: %d bytes\n", sizeof(MdaData));
    printf("MdaDataStruct is at: %p\n", temp);
    printf("Number of points is: %d\n", numPts);
    printf("Number of L values is: %d\n", numL);
    printf("Data size is: %d bytes\n", sizeof(float)*numPts);
    printf("Data is at: %p\n", temp->data);
    printf("Distribution size is: %d bytes\n", sizeof(float)*numPts*numL);
    printf("Distributions are at: %p\n", temp->dists);
    printf("Residual size is: %d bytes\n", sizeof(float)*numPts);
    printf("Residuals are at: %p\n", temp->resids);
#endif
    return ((void*)temp);
}

void freeMdaStruct(void* strPtr)
{
    MdaData* temp = (MdaData*)strPtr;
#ifdef MAKE_AND_FREE_OUTPUT
    printf("\n\nFreeing the residuals at %p\n", temp->resids);
#endif
    free((temp->resids));
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the distributions at %p\n", temp->dists);
#endif
    free((temp->dists));
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the data at %p\n", temp->data);
#endif
    free((temp->data));
#ifdef MAKE_AND_FREE_OUTPUT
    printf("Freeing the structure at %p\n", temp);
#endif
    free(temp);
}

void setMdaData(void* strPtr, float* divData)
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
        printf("At cell %d value: %f | Orig: %f\n", i, data->data[i], divData[i]);
#endif
    }
}

void setMdaDist(void* strPtr, int distIndex, float* dist)
{
    MdaData* data = (MdaData*)strPtr;
    if (distIndex >= data->numLs)
    { return; }

#ifdef ASSIGN_DIST_OUTPUT
    printf("struct ptr is: %p\n", data);
    printf("dist ptr is: %p\n", data->dists);
    printf("index is: %d\n",distIndex);
#endif
    
    int offset = distIndex*data->numPts;
#ifdef ASSIGN_DIST_OUTPUT
    printf("Offset is: %d\n", offset);
#endif
    
    float* temp = &(data->dists[offset]);
    
    for (int i=0; i<data->numPts; ++i)
    {
        //data->dists[offset+i] = dist[i];
        temp[i] = dist[i];
#ifdef ASSIGN_DIST_OUTPUT
        printf("At cell %d value: %f  | Orig: %f\n", i, data->dists[i], dist[i]);
#endif
    }
}

float calculateChi(void* strPtr, float* params)
{
    MdaData* data = (MdaData*)strPtr;
    int numPts = data->numPts;
    int numL = data->numLs;
    float* exp = data->data;
    float* dists = data->dists;
    float* res = data->resids;
    float pVal = params[0];
    float* dist = dists;
    
    //first use only the first distribution to load the residuals
    for (int j=0; j<numPts; ++j)
    {
        res[j] = (exp[j]-(pVal*dist[j]));
    }

    //iterate across the distributions
    for(int i=1; i<numL; ++i)
    {
        float pVal = params[i];
        float* dist = &(dists[i*numPts]);
        for (int j=0; j<numPts; ++j)
        {
            res[j] -= (pVal*dist[j]);
        }
    }
    
    float chi = 0.0f;
    
    //now accumulate the residuals squared and return them
    for(int i=0; i<numPts; ++i)
    {
        chi += (res[i]*res[i]);
    }
    
    return chi;
}


//calculate the log liklihood
float calculateLnLiklihood(void* strPtr, float* params)
{
    MdaData* data = (MdaData*)strPtr;
    int numPts = data->numPts;
    int numL = data->numLs;
    float* exp = data->data;
    float* dists = data->dists;
    float* res = data->resids;
    float pVal = params[0];
    float* dist = dists;
    
    //subtract the first distribution to load the residuals
    for (int j=0; j<numPts; ++j)
    {
        res[j] = (exp[j]-(pVal*dist[j]));
    }

    //iterate across the distributions
    for(int i=1; i<numL; ++i)
    {
        float pVal = params[i];
        float* dist = &(dists[i*numPts]);
        for (int j=0; j<numPts; ++j)
        {
            res[j] -= (pVal*dist[j]);
        }
    }
    
    float chi = 0.0f;
    
    //now accumulate the residuals squared and return them
    for(int i=0; i<numPts; ++i)
    {
        chi += (res[i]*res[i]);
    }
    
    return (chi/(-2.0));
}

//calculate the log liklihood
float calculateLnLiklihoodResids(void* strPtr, float* params, float* residArray)
{
    MdaData* data = (MdaData*)strPtr;
    int numPts = data->numPts;
    int numL = data->numLs;
    float* exp = data->data;
    float* dists = data->dists;
    float pVal = params[0];
    float* dist = dists;
    
    //subtract the first distribution to load the residuals
    for (int j=0; j<numPts; ++j)
    {
        residArray[j] = (exp[j]-(pVal*dist[j]));
    }

    //iterate across the distributions
    for(int i=1; i<numL; ++i)
    {
        float pVal = params[i];
        float* dist = &(dists[i*numPts]);
        for (int j=0; j<numPts; ++j)
        {
            residArray[j] -= (pVal*dist[j]);
        }
    }
    
    float chi = 0.0f;
    
    //now accumulate the residuals squared and return them
    for(int i=0; i<numPts; ++i)
    {
        chi += (residArray[i]*residArray[i]);
    }
    
    return (chi/(-2.0));
}
