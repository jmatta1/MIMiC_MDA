#include"ChiSquare.h"
#include<stdio.h>
#include<time.h>

int main(int argc, char* argv[])
{
    clock_t t1, t2;
    t1=clock();
    int diff;
    MdaData* dataPtr = makeMdaStruct(10, 3);
    freeMdaStruct(dataPtr);
    //test the allocate function
    t1 = clock();
    dataPtr = makeMdaStruct(10, 3);
    t2 = clock();
    diff = (((t2-t1)*1000000)/CLOCKS_PER_SEC);
    printf("make struct, took: %d microseconds\n", diff);
    
    //make some data
    float testArr[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    
    //test the set data function
    t1 = clock();
    setMdaData(dataPtr, testArr);
    t2 = clock();
    diff = (((t2-t1)*1000000)/CLOCKS_PER_SEC);
    printf("set data, took: %d microseconds\n", diff);
    
    //test the set dist function
    t1 = clock();
    setMdaDist(dataPtr, 0, testArr);
    setMdaDist(dataPtr, 1, testArr);
    setMdaDist(dataPtr, 2, testArr);
    t2 = clock();
    diff = (((t2-t1)*1000000)/CLOCKS_PER_SEC);
    printf("set dist x3, took: %d microseconds\n", diff);
    
    //high res timers for the looping
    struct timespec t3;
    struct timespec t4;
    
    //test the calculate chi function
    float paramArray[3] = {0.2, 0.2, 0.2};
    clock_gettime(CLOCK_MONOTONIC , &t3);
    float chi = calculateChi(dataPtr, paramArray);
    for( int i=0; i< 30; ++i)
    {
        calculateChi(dataPtr, paramArray);
    }
    clock_gettime(CLOCK_MONOTONIC , &t4);
    diff = (t4.tv_nsec-t3.tv_nsec);
    printf("calc chi, took: %d nanoseconds\n", diff);
    printf("  the chi was: %f\n", chi);
    
    //test the deallocate function
    t1 = clock();
    freeMdaStruct(dataPtr);
    t2 = clock();
    diff = (((t2-t1)*1000000)/CLOCKS_PER_SEC);
    printf("free struct, took: %d microseconds\n", diff);
    
    //exit
    return 0;
}
