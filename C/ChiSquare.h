#ifndef MDA_MCMC_C_CHI_SQUARE_H
#define MDA_MCMC_C_CHI_SQUARE_H

// define the structure that holds the various pointers and values used to
// store the data (divided by the errors) and the interpolated distribution
// points divided by the errors. In this fashion all the division is taken
// care of ahead of time and instead all that needs to happen is the function
// be summed from the predivided interpolated points times the appropriate sum
// rule fraction, then those points are subtracted from the data, the residuals
// are summed and finally the squared residuals are summed, this sum is then
// divided by two and made negative to yield the log likelihood.
// once this structure is created, the distributions loaded and the data loaded
// it is not modified by the likelihood function calls
typedef struct MdaDataStruct
{
    //stores the number of data points
    int numPts;
    //stores the number of L values
    int numLs;
    //stores the predivided data points
    //contains numPoints cells
    float* data;
    //stores the set of interpolated distributions
    //contains numPoints*numLVals cells
    float* dists;
    //stores the list of residuals needed on a temporary basis
    float* resids;
} MdaData;

//allocate the structure
__attribute__((malloc)) void* makeMdaStruct(int numPts, int numL);

//deallocate the structure
void freeMdaStruct(void* strPtr);

//assign the data array
void setMdaData(void* strPtr, float* divData);

//assign one of the distributions
void setMdaDist(void* strPtr, int distIndex, float* dist);

//calculate the chi^2
float calculateChi(void* strPtr, float* params);

#endif  // MDA_MCMC_C_CHI_SQUARE_H
