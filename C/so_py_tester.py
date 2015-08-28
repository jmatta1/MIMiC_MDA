from ctypes import *
import numpy as np

def main():
    # load the library
    chi_sq_lib = cdll.LoadLibrary("./libChiSq.so")
    # extract the function to generate the struct
    makeMdaStruct = chi_sq_lib.makeMdaStruct
    # extract the function to free the struct
    freeMdaStruct = chi_sq_lib.freeMdaStruct
    # tell python the function will return a void pointer
    makeMdaStruct.restype = c_void_p
    # make the strcuture to hold the MdaData with 50 data points and 8 L dists
    mdaData = makeMdaStruct(50,8)
    # print the pointer to the structure
    print hex(mdaData)
    #free the structure
    freeMdaStruct(mdaData)



if __name__ == "__main__":
    main()
