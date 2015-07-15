import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import ctypes
import time


def rmod(x, y):
    result = x - int(x/y)*y
    return result if result >= 0 else result + y


start=time.time()
LIB = '/Users/martin/Documents/chaos/ccodefromsage/invariantdynamics/ergodic.dylib'

grid=410
imatrix = np.ones((grid,grid))#np.loadtxt("outputs/imatrix.txt")
imatrix = imatrix.astype(long)
#p = np.loadtxt("outputs/parameters.txt")
xmin=0.
ymin=0.
xmax=2*np.pi;
ymax=2*np.pi
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

kgrid = 25;
totaltime = 1000;
ERG = ctypes.CDLL(LIB)
print "started c"
"""ERG.partition(ctypes.c_int32(int(p[0])),
        ctypes.c_int32(int(p[0])),
        ctypes.c_int32(totaltime),
        ctypes.c_double(p[1]),
        ctypes.c_double(p[2]),
        ctypes.c_double(p[3]),
        ctypes.c_double(p[4]),
        ctypes.c_int32(kgrid),
        imatrix.ctypes.data_as(ctypes.c_void_p))"""
ERG.partition(ctypes.c_int32(grid),
        ctypes.c_int32(grid),
        ctypes.c_int32(totaltime),
        ctypes.c_double(xmin),
        ctypes.c_double(ymin),
        ctypes.c_double(deltax),
        ctypes.c_double(deltay),
        ctypes.c_int32(kgrid),
        imatrix.ctypes.data_as(ctypes.c_void_p))
print "finished c"
w = np.vectorize(lambda x: x)(imatrix)
imgplot = plt.imshow(w,interpolation='none' if grid > 200 else 'nearest',cmap=cm.flag,extent=[xmin,xmin+grid*deltax,ymin,ymin+grid*deltay])
plt.savefig("outputs/ergodic_result.ps")
plt.clf()
end = time.time()
