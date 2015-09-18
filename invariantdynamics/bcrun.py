import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
CBIRK=CDLL(liblist[13].strip())
libfile.close()

#Parameters


#Number of Functions Converging
fnum = 1

#Curve Parameters
ax = 0
bx = np.pi
ay = 0
by = 0
#Weighted?
weighted=0
#Iteration Number
totaltime = 10000   # Map iterates
numpoints = 1000 #sampling of curve
#Save Parameters
m = np.zeros((numpoints,2),dtype="float64")
#print "started birk_curve.c"
CBIRK.curve_part(c_int32(weighted),
        c_double(ax),
        c_double(ay),
        c_double(bx),
        c_double(by),
        c_int32(numpoints),
        c_int32(totaltime),
        m.ctypes.data_as(c_void_p))
    #print "finished birk_curve.c"
m = np.transpose(m)
plt.plot(m[0]*np.pi,m[1],'-')
plt.xlabel('y')
plt.ylabel('phi(y)')
plt.savefig('outputs/birk_curve_result.pdf')
plt.clf()
