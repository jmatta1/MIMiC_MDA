from ctypes import *
import numpy as np
import time

def main():
    print "starting python testing"
    # load the library
    cs_lib = cdll.LoadLibrary("./libChiSq.so")
    # tell python the function will return a void pointer
    cs_lib.makeMdaStruct.restype = c_void_p
    # make the strcuture to hold the MdaData with 50 data points and 8 L dists
    t1 = time.time()
    mdaData = cs_lib.makeMdaStruct(10,3)
    t2 = time.time()
    print "make struct, took", (t2-t1), "seconds"
    # make an array with the test data
    divDat = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], np.float32)
    #run the function to load the array into the struct
    t1 = time.time()
    cs_lib.setData(mdaData, divDat.ctypes.data)
    t2 = time.time()
    print "set data, took", (t2-t1), "seconds"
    #free the structure
    t1 = time.time()
    cs_lib.freeMdaStruct(mdaData)
    t2 = time.time()
    print "free struct, took", (t2-t1), "seconds"


if __name__ == "__main__":
    main()
