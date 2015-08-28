#include"ChiSquare.h"
#include<stdio.h>
#include<time.h>

int main(int argc, char* argv[])
{
    clock_t t1, t2;
    t1=clock();
    float diff;
    printf("starting C testing\n");
    //test the allocate function
    t1 = clock();
    MdaData* dataPtr = makeMdaStruct(10, 3);
    t2 = clock();
    diff = (((float)(t2-t1))/CLOCKS_PER_SEC);
    printf("make struct, took: %f\n", diff);
    
    //make some data
    float testArr[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    
    //test the set data function
    t1 = clock();
    setData(dataPtr, testArr);
    t2 = clock();
    diff = (((float)(t2-t1))/CLOCKS_PER_SEC);
    printf("set data, took: %f\n", diff);
    
    //test the deallocate function
    t1 = clock();
    freeMdaStruct(dataPtr);
    t2 = clock();
    diff = (((float)(t2-t1))/CLOCKS_PER_SEC);
    printf("del struct, took: %f\n", diff);
    
    //exit
    return 0;
}
