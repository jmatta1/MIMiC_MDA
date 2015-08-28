#include"ChiSquare.h"
#include<stdio.h>
#include<time.h>

int main(int argc, char* argv[])
{
    clock_t t1, t2;
    t1=clock();
    long long int diff;
    printf("starting C testing\n");
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
    
    //test the deallocate function
    t1 = clock();
    freeMdaStruct(dataPtr);
    t2 = clock();
    diff = (((t2-t1)*1000000)/CLOCKS_PER_SEC);
    printf("free struct, took: %d microseconds\n", diff);
    
    //exit
    return 0;
}
