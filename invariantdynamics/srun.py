import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import cPickle as pickle
import time

start = time.time()

#Required Libraries
libfile = open("../libraries.txt")
liblist = libfile.readlines()
STICKY=CDLL(liblist[5].strip())
libfile.close()


#Parameters

#Number of Windows
windows = 400;

#Closeness Requirement
epsilon = 0.5;

#Initial Conditions
ix = 5.14;
iy = 2.14;

rrs = np.zeros(windows, dtype="float")

print "started stickiness.c"
STICKY.stickiness(c_int32(windows),
        c_double(ix),
        c_double(iy),
        c_double(epsilon),
        rrs.ctypes.data_as(c_void_p))
print "finished stickiness.c"

print rrs

print "Plot Recurrence Rate"
plt.plot(rrs)
plt.ylabel("Recurrence Rate")
plt.xlabel("Window")
plt.savefig("outputs/rr_result.ps")
plt.clf()

end=time.time()
print(end-start)
