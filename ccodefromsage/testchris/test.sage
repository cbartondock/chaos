import numpy
import matplotlib.pylab as plt
import ctypes
import time
start=time.time()
LIB = '/Users/chrisdock/Documents/chaos/ccodefromsage/testchris/test.dylib'

param1 = 1.4
param2 = .3
maxiter = 10
numsamples = 200

"""
ymin=0
ymax=2*N(pi)
xmin=0
xmax=2*N(pi)
"""
ymin = -2.5
ymax = 2.5 
xmin = -2.5
xmax = 2.5

grid=1000

deltax = (xmax-xmin)/grid
deltay = (ymax-ymin)/grid

mapnum = 1


X = ctypes.CDLL(LIB)
m = matrix.ones(grid,grid)
m = m.numpy('uint8')
print "started c"
X.forward(ctypes.c_ubyte(mapnum), 
       ctypes.c_double(param1),
       ctypes.c_double(param2),
	   ctypes.c_int32(maxiter),
	   ctypes.c_int32(numsamples),
       ctypes.c_int32(grid),
       ctypes.c_int32(grid),
	   ctypes.c_double(xmin),
	   ctypes.c_double(ymin),
	   ctypes.c_double(deltax),
       ctypes.c_double(deltay),
	   m.ctypes.data_as(ctypes.c_void_p))
print "finished c"
x = numpy.vectorize(lambda x: x)(m)
imgplot = plt.imshow(x,vmin=0, vmax=1,interpolation='nearest')
plt.savefig("algc_result.png")
end=time.time()
print "total time: "+str(end-start)


