import numpy
import ctypes
import os

LIB = '/Users/chrisdock/Documents/chaos/ccodefromsage/testvincent/test.dylib'

n = 50      # 100 x 100 resulution
S = 200      # samples per pixel
I = 10       # iterations of map

rho = 0.3
c1  = 1.4

ll_x = -2.5  # lower left X
ll_y = -2.5  # lower left Y

dx   = 5.0/n # pixel size


X = ctypes.CDLL(LIB)
m = matrix.ones(n,n)
m = m.numpy('int64')

X.forward( ctypes.c_double(c1),
           ctypes.c_double(rho),
	   ctypes.c_int32(I),
	   ctypes.c_int32(S),
	   ctypes.c_double(ll_x),
	   ctypes.c_double(ll_y),
	   ctypes.c_double(dx),
	   ctypes.c_int32(n),
	   ctypes.c_int32(n),
	   m.ctypes.data_as(ctypes.c_void_p))

x = numpy.vectorize(lambda x: x)(m)
matrix_plot(x,cmap='gray_r')



