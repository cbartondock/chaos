import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
STICKY=CDLL(liblist[5].strip())
libfile.close()


#Parameters
windowwidth=1000;
epsilon=.5;
ix = 5.14;
iy=2.14;

m = np.zeros((windowwidth,windowwidth),dtype="uint32")
xlist = np.zeros(windowwidth, dtype="float")
ylist = np.zeros(windowwidth, dtype="float")


STICKY.rp_window(c_int32(windowwidth),
        c_double(ix),
        c_double(iy),
        c_double(epsilon),
        m.ctypes.data_as(c_void_p),
        xlist.ctypes.data_as(c_void_p),
        ylist.ctypes.data_as(c_void_p))
if windowwidth <= 50:
    print(m)
plt.imshow(m,vmin=0,interpolation='nearest',cmap=cm.Reds,extent=[0,windowwidth,0,windowwidth])
plt.colorbar()
plt.savefig("outputs/rp_result.ps")
plt.clf()
end=time.time()
print(end-start)
