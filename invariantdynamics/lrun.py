import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
LYAP=CDLL(liblist[10].strip())
libfile.close()

#Sample points
sx, sy = 1000, 1000

points = [[2*np.pi*(x+.5)/float(sx), 2*np.pi*(y+.5)/float(sy),0.] for x in range(0,sx) for y in range(0,sy)]
points = np.array(points)
pointspt = (POINTER(c_double)*len(points))(*[row.ctypes.data_as(POINTER(c_double)) for row in points])

totaltime = 1000 # Map iterates


print "started extended_smap.c"
LYAP.exponents(pointspt,
        c_int32(len(points)),
        c_int32(totaltime))
print "finished extended_smap.c"

m = [[0 for i in range(0,sx)] for j in range(sy)]
print
for point in points:
    m[int(point[1]*sy/(2*np.pi))][int(point[0]*sx/(2*np.pi))] = point[2]

plt.imshow(m)
plt.colorbar()
plt.savefig("outputs/lyapunovplot.png")
plt.clf()
end = time.time()
print(end-start)
