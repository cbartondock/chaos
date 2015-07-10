import numpy
import matplotlib.pylab as plt
import matplotlib.cm as cm
import ctypes
import time
start=time.time()
LIB = '/Users/chrisdock/Documents/chaos/ccodefromsage/testfast/test.dylib'

maxiter = 500
numsamples = 100
grid= 128
print "grid: "+str(grid)+", numsamples: "+str(numsamples)+", maxiter: "+str(maxiter)

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

deltax = (xmax-xmin)/grid
deltay = (ymax-ymin)/grid


X = ctypes.CDLL(LIB)
m = matrix.ones(grid,grid)
m = m.numpy('uint8')
print "started c"
X.forward(ctypes.c_int32(maxiter),
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
imgplot = plt.imshow(x,vmin=0, vmax=1,interpolation='nearest',cmap=cm.winter,extent=[xmin,xmax,ymin,ymax])
plt.savefig("algc_result.png")
end=time.time()
print "total time: "+str(end-start)


