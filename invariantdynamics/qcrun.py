
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
CQUAS=CDLL(liblist[7].strip())
libfile.close()

#Parameters


#Number of Functions Converging
fnum = 1

#Curve Parameters
ax = 0
bx = 6.28
ay = 6.28
by = 6.28


#Iteration Number
totaltime = 10000   # Map iterates
numpoints = 100 #sampling of curve
#Save Parameters

m = np.zeros((numpoints,2),dtype="float64")
print "started quasiperiodicity.c"
CQUAS.curve_convergence(c_double(ax),
                 c_double(ay),
                 c_double(bx),
                 c_double(by),
                 c_int32(numpoints),
                 c_int32(totaltime),
		         c_int32(fnum),
                 m.ctypes.data_as(c_void_p))
print "finished quasiperiodicity.c"
m = np.transpose(m)
print m[0]
print m[1]

plt.plot(m[0],m[1],'-')




plt.xlabel('parameter')
plt.ylabel('#zeros')
plt.savefig('outputs/quasi_curve_result.pdf')

"""
#Plotting Convergence Spatially
w = np.vectorize(lambda x: x)(m)
print(m)
plt.imshow(w,vmin=0,interpolation='nearest',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
cbar = plt.colorbar()
#plt.title("Rates of Birkhoff Convergence in the Standard Map (N={0})".format(totaltime))
cbar.set_label("#zeros")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("outputs/convergence_result.ps")
plt.clf()
end=time.time()
print(end-start)
"""
