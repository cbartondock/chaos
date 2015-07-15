import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

start = time.time()
QUAS=CDLL("/Users/chrisdock/Documents/chaos/ccodefromsage/invariantdynamics/quasiperiodicity.dylib")

grid=500
m = np.ones((grid,grid))#np.loadtxt("outputs/imatrix.txt")
m = m.astype('uint8')
#p = np.loadtxt("outputs/parameters.txt")
xmin=0.
ymin=0.
xmax=2*np.pi;
ymax=2*np.pi
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid
totaltime = 10000

QUAS.convergence(c_int32(grid),
                 c_int32(grid),
                 c_int32(totaltime),
                 c_double(xmin),
                 c_double(ymin),
                 c_double(deltax),
                 c_double(deltay),
                 m.ctypes.data_as(c_void_p))
w = np.vectorize(lambda x: x)(m)
print(m)
plt.imshow(w,interpolation='none' if grid > 200 else 'nearest',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
plt.colorbar()
plt.savefig("outputs/convergence_result.ps")

end=time.time()
print(end-start)
