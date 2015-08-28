#include"ChiSquare.h"
#include<stdio.h>

int main(int argc, char* argv[])
{
    MdaData* dataPtr = makeMdaStruct(50, 8);
    printf("%p\n",dataPtr);
    freeMdaStruct(dataPtr);
    return 0;
}
