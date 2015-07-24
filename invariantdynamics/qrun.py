import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
QUAS=CDLL(liblist[4].strip())
libfile.close()

'''
Parameters
'''
grid=500           # Grid size nxn
fnum = 1           # Number of functions
xmin=0.            #
ymin=0.            # Dimensions
xmax=2*np.pi       #
ymax=2*np.pi       #
totaltime = 10000   # Map iterates
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

np.savetxt("outputs/qparameters.txt", np.array([int(grid), float(xmin),float(ymin),float(deltax),float(deltay)]))

m = np.ones((grid,grid),dtype="uint8")

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
entries = m.flatten()
plt.hist(entries,np.arange(min(entries),max(entries),1)-1.5)
plt.savefig("outputs/convergence_histogram.ps")

outfile = open("outputs/conv_matrix.p", "wb")
pickle.dump(m,outfile)
outfile.close()
end=time.time()
print(end-start)