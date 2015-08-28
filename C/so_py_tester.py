from ctypes import cdll

chi_sq_lib = cdll.LoadLibrary("./C/libChiSq.so")

hello = chi_sq_lib.hello

print hello("world")

from ctypes import c_char_p

hello.restype = c_char_p

print hello("world")
