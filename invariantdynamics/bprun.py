import numpy as np
import os
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time

def droplist(l):
    l.sort()
    mydict = {}
    stepper = 0
    for item in l:
        if item not in mydict:
            mydict[item] = stepper
            stepper += 1
    return mydict
start=time.time()

cmapname = "jet"

#Required Libraries
libraries = open('../libraries.txt').readlines()
ERG = CDLL(libraries[3].strip())

#Parameters (running on the Standard Map)
grid=1500

#Region Parameters
xmin=1.
ymin=0
xmax=np.pi
ymax=np.pi-1.

deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

#Partition Sensitivity
kgrid = 50

#Iterations
totaltime = 1000

#Weighted?
weighted=0

m = np.ones((grid,grid),dtype="uint64")

print "started birkhoff_partition.c"
ERG.partition(c_int32(weighted),
        c_int32(grid),
        c_int32(grid),
        c_int32(totaltime),
        c_double(xmin),
        c_double(ymin),
        c_double(deltax),
        c_double(deltay),
        c_int32(kgrid),
        m.ctypes.data_as(c_void_p))
print "finished birkhoff_partition.c"
colordict = droplist([entry for row in m for entry in row])
m = [[colordict[entry] for entry in row] for row in m]
w = np.vectorize(lambda x: x)(m)
#Plotting
plt.imshow(m,interpolation='none' if grid > 200 else 'nearest',cmap=cm.get_cmap(cmapname),extent=[xmin,xmin+grid*deltax,ymin+grid*deltay,ymin])
#plt.title("Invariant Curves of the Standard Map")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("outputs/birkhoff_result_t{0}_g{1}_xs{2}_ys{3}_xb{4}_yb{5}_{6}_2.pdf".format(totaltime,grid,round(xmin,2),round(ymin,2),round(xmax,2),round(ymax,2), "weighted" if weighted else "unweighted"))
plt.clf()
end = time.time()
print end-start

weighted = 1

m = np.ones((grid,grid),dtype="uint64")

print "started birkhoff_partition.c"
ERG.partition(c_int32(weighted),
        c_int32(grid),
        c_int32(grid),
        c_int32(totaltime),
        c_long_double(xmin),
        c_long_double(ymin),
        c_long_double(deltax),
        c_long_double(deltay),
        c_int32(kgrid),
        m.ctypes.data_as(c_void_p))
print "finished birkhoff_partition.c"
colordict = droplist([entry for row in m for entry in row])
m = [[colordict[entry] for entry in row] for row in m]
w = np.vectorize(lambda x: x)(m)
#Plotting
plt.imshow(m,interpolation='none' if grid > 200 else 'nearest',cmap=cm.get_cmap(cmapname),extent=[xmin,xmin+grid*deltax,ymin+grid*deltay,ymin])
#plt.title("Invariant Curves of the Standard Map")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("outputs/birkhoff_result_t{0}_g{1}_xs{2}_ys{3}_xb{4}_yb{5}_{6}_2.pdf".format(totaltime,grid,round(xmin,2),round(ymin,2),round(xmax,2),round(ymax,2), "weighted" if weighted else "unweighted"))
plt.clf()
end = time.time()
print end-start
