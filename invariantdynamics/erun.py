import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
SEARCH=CDLL(liblist[8].strip())
libfile.close()

#Sample points
sx, sy, sz = 100, 100, 100
points = [[x/float(sx),y/float(sy),z/float(sz),0.] for x in range(0, sx) for y in range(0,sy) for z in range(0,sz)]
#points = [[0.499 +.002*x/float(sx), 0.,0.,0.] for x in range(0,sx)]
points = np.array(points)
pointspt = (POINTER(c_double)*len(points))(*[row.ctypes.data_as(POINTER(c_double)) for row in points])

totaltime = 100  # Map iterates


print "started extended_smap.c"
SEARCH.quasi_search(pointspt,
        c_int32(len(points)),
        c_int32(totaltime))
print "finished extended_smap.c"


#Plotting Convergence Rates Histogram
points = [[pointspt[i][k] for k in range(0,4)] for i in range(0,len(points))]
points = np.array([point for point in points if point[3]!=float('Inf')])
convdata = np.transpose(points)
plt.hist(convdata[3], np.arange(min(convdata[3]),max(convdata[3]),.1),alpha=0.6)



interesting_points =[]
for convpoint in points:
    if convpoint[3]>20:
        interesting_points.append(convpoint)
print interesting_points
plt.savefig("outputs/esmap_crates_hist.png")
""""m = [[0 for i in range(0,sx)] for j in range(sy)]
for point in points:
    m[int(point[0]*1000)][int(point[1]*1000)] = point[3]

plt.imshow(m)
plt.colorbar()
plt.show()
"""
"""
fig = plt.figure()
ax = Axes3D(fig)
p = ax.scatter(xs=convdata[0],ys=convdata[1],zs=convdata[2],s=36,c=convdata[3])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
fig.colorbar(p)
plt.savefig("outputs/esmap_crates_plot.png")
"""
end = time.time()
print(end-start)
