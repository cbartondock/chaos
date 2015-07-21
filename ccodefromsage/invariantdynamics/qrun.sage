import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

start = time.time()

libfile = open("libraries.txt")
liblist = libfile.readlines()
QUAS=CDLL(liblist[4].strip())
libfile.close()

'''
Parameters
'''
grid=1000           # Grid size nxn
fnum = 3           # Number of functions
xmin=0.            #
ymin=0.            # Dimensions
xmax=2*np.pi       #
ymax=2*np.pi       #
totaltime = 10000   # Map iterates

m = np.ones((grid,grid))#np.loadtxt("outputs/imatrix.txt")
m = m.astype('uint8')
#p = np.loadtxt("outputs/parameters.txt")

deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

QUAS.convergence(c_int32(grid),
                 c_int32(grid),
                 c_int32(totaltime),
                 c_double(xmin),
                 c_double(ymin),
                 c_double(deltax),
                 c_double(deltay),
		 c_int32(fnum),
                 m.ctypes.data_as(c_void_p))
w = np.vectorize(lambda x: x)(m)
print(m)
plt.imshow(w,vmin=0,interpolation='nearest',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
plt.colorbar()
plt.savefig("outputs/convergence_result.ps")
plt.clf()
plt.hist(m.flatten(), 50, histtype='stepfilled')
plt.savefig("outputs/convergence_histogram.ps")
end=time.time()
print(end-start)
