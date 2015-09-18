import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
SEARCH=CDLL(liblist[9].strip())
libfile.close()

#4.514155 (100)
#4.54614, .00002 (10000)
#Sample points
sx, sy = 100, 100
w=.000000000001
yval = 1#1.7499999

#points = [[2*np.pi*(x+.5)/float(sx), 2*np.pi*(y+.5)/float(sy),0.] for x in range(0,sx) for y in range(0,sy)]
points = [[3.14 + w*((x+.5)/float(sx)), yval + w*((y+.5)/float(sy)),0.] for x in range(0, sx) for y in range(0,sy)]
points = np.array(points)
pointspt = (POINTER(c_double)*len(points))(*[row.ctypes.data_as(POINTER(c_double)) for row in points])

totaltime = 10000 # Map iterates


print "started extended_smap.c"
SEARCH.quasi_search(pointspt,
        c_int32(len(points)),
        c_int32(totaltime),
        c_int32(0))
print "finished extended_smap.c"
m = [[0 for i in range(0,sx)] for j in range(sy)]
print
for point in points:
    m[int((point[1]-yval)*sy/w)][int((point[0]-3.14)*sx/w)] = point[2]
    #m[int(point[1]*sy/(2*np.pi))][int(point[0]*sx/(2*np.pi))] = point[2]

plt.imshow(m,extent=[3.14,3.14+.0000000001,1.7499999,1.7499999+.0000000001])
plt.colorbar()
#plt.show()
plt.savefig("outputs/cm_weightless_partition2.png")
plt.clf()

testpoint = m[int(sy/2.)][int(sx/2.)]
m2 = [[1 if m[i][j] >= testpoint else 0 for j in range(sx)] for i in range(sy)]
plt.imshow(m2,vmin=0,vmax=1,extent=[3.14,3.14+w,yval,yval+w])
#plt.show()
plt.savefig("outputs/cm_weightless_divide2.png")



end = time.time()
print(end-start)
