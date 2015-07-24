import numpy as np
import os
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

start=time.time()
libraries = open('../libraries.txt').readlines()
ERG = CDLL(libraries[3].strip())

grid=150
imatrix = np.ones((grid,grid),dtype="uint64")

xmin=0.
ymin=0.
xmax=2*np.pi;
ymax=2*np.pi
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

kgrid = 25;
totaltime = 1000;
print "started c"
ERG.partition(c_int32(grid),
        c_int32(grid),
        c_int32(totaltime),
        c_double(xmin),
        c_double(ymin),
        c_double(deltax),
        c_double(deltay),
        c_int32(kgrid),
        imatrix.ctypes.data_as(c_void_p))
print "finished c"
w = np.vectorize(lambda x: x)(imatrix)
imgplot = plt.imshow(w,interpolation='none' if grid > 200 else 'nearest',cmap=cm.flag,extent=[xmin,xmin+grid*deltax,ymin,ymin+grid*deltay])
plt.savefig("outputs/ergodic_result.ps")
plt.clf()
end = time.time()