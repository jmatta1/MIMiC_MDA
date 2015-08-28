from ctypes import *

chi_sq_lib = cdll.LoadLibrary("./libChiSq.so")

makeMdaStruct = chi_sq_lib.makeMdaStruct

makeMdaStruct.restype = c_void_p

mdaData = makeMdaStruct(50,8)

print mdaData
