from ctypes import *
import numpy as np
import time

def main():
    
    # load the library
    t1 = time.time()
    cs_lib = cdll.LoadLibrary("./libChiSq.so")
    t2 = time.time()
    print "load library, took", int((t2-t1)*1000000), "microseconds"
    
    # tell python the function will return a void pointer
    cs_lib.makeMdaStruct.restype = c_void_p
    cs_lib.calculateChi.restype = c_float
    cs_lib.calculateLnLiklihood.restype = c_float
    
    # make the strcuture to hold the MdaData with 50 data points and 8 L dists
    t1 = time.time()
    mdaData = cs_lib.makeMdaStruct(10,3)
    t2 = time.time()
    print "make struct, took", int((t2-t1)*1000000), "microseconds"
    
    # make an array with the test data
    divDat = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], np.float32)
    
    # run the function to load the array into the struct
    t1 = time.time()
    cs_lib.setMdaData(mdaData, divDat.ctypes.data)
    t2 = time.time()
    print "set data, took", int((t2-t1)*1000000), "microseconds"
    
    # run the function to load the array into the struct
    t1 = time.time()
    cs_lib.setMdaDist(mdaData, 0, divDat.ctypes.data)
    cs_lib.setMdaDist(mdaData, 1, divDat.ctypes.data)
    cs_lib.setMdaDist(mdaData, 2, divDat.ctypes.data)
    t2 = time.time()
    print "set dist x3, took", int((t2-t1)*1000000), "microseconds"
    
    # make an array with the test data
    params = np.array([0.2, 0.2, 0.2], np.float32)
    
    # run the function to calculate a chi^2
    t1 = time.time()
    chi = cs_lib.calculateChi(mdaData, params.ctypes.data)
    for _ in range(30):
        cs_lib.calculateChi(mdaData, params.ctypes.data)
    t2 = time.time()
    print "calc chi, took", int((t2-t1)*1000000), "microseconds"
    print "  the chi was:",chi
    
    # run the function to calculate a chi^2
    t1 = time.time()
    lnlik = cs_lib.calculateLnLiklihood(mdaData, params.ctypes.data)
    for _ in range(30):
        cs_lib.calculateLnLiklihood(mdaData, params.ctypes.data)
    t2 = time.time()
    print "calc log liklihood, took", int((t2-t1)*1000000), "microseconds"
    print "  the log liklihood was:",lnlik
    
    # free the structure
    t1 = time.time()
    cs_lib.freeMdaStruct(mdaData)
    t2 = time.time()
    print "free struct, took", int((t2-t1)*1000000), "microseconds"
    
    # as a test, calculate the chi^2 in the pythonic way
    t1 = time.time()
    residuals = (divDat - divDat*params[0])
    for param in params[1:]:
        residuals -= param*divDat
    lnlik = np.sum(np.power(residuals,2.0))/(-2.0)
    for _ in range(30):
        residuals = (divDat - divDat*params[0])
        for param in params[1:]:
            residuals -= param*divDat
        temp = np.sum(np.power(residuals,2.0))/(-2.0)
    t2 = time.time()
    print "python calc log liklihood, took", int((t2-t1)*1000000), "microseconds"
    print "  the log liklihood was:",lnlik


if __name__ == "__main__":
    main()
