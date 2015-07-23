import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import cPickle as pickle
import time

start = time.time()

libfile = open("libraries.txt")
liblist = libfile.readlines()
STICKY=CDLL(liblist[5].strip())
libfile.close()


#Parameters
windows = 400;
epsilon = 0.5;
ix = 5.14;
iy = 2.14;

rrs = np.zeros(windows)
rrs = rrs.astype('float')
STICKY.stickiness(c_int32(windows),
        c_double(ix),
        c_double(iy),
        c_double(epsilon),
        rrs.ctypes.data_as(c_void_p))

print rrs
plt.plot(rrs)
plt.ylabel("Recurrence Rate")
plt.xlabel("Window")
plt.savefig("outputs/rr_result.ps")
plt.clf()
end=time.time()
print(end-start)
