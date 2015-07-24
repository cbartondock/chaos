import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from ctypes import *
import time
start=time.time()
libraries = open("../libraries.txt").readlines()
FAC = CDLL(libraries[6].strip())
doublify = lambda x: '0'*(64-len(bin(x)[2:]))+bin(x)[2:]
maxiter = 30
numsamples = 30
grid= 100
mul=2
print "grid: "+str(grid)+", numsamples: "+str(numsamples)+", maxiter: "+str(maxiter)
"""
ymin=0
ymax=2*N(pi)
xmin=0
xmax=2*N(pi)
"""
ymin = -2.
ymax = 2.
xmin = -2.
xmax = 2.
"""
ymin=-2.8
ymax=1.8
xmin=-1.8
xmax=2.8
"""
deltax = (xmax-xmin)/grid
deltay = (ymax-ymin)/grid

m = np.ones((grid,grid),dtype="uint64")
m = np.multiply(2, m)

h = np.zeros((mul*grid,mul*grid),dtype="uint64");
b = [[0 for i in range(0,grid)] for j in range(0,grid)]
for i in range(0,grid):
    for j in range(0,grid):
        if i==0 or j==0 or i==grid-1 or j==grid-1:
            b[i][j]=1
x3 = np.vectorize(lambda x:int(x)/1.)(b)
x3 = np.ma.masked_where(x3 <0.9,x3)
print "started c"
FAC.generate_images(c_int32(mul),
       c_int32(maxiter),
       c_int32(maxiter),
	   c_int32(numsamples),
       c_int32(grid),
       c_int32(grid),
	   c_double(xmin),
	   c_double(ymin),
	   c_double(deltax),
       c_double(deltay),
	   m.ctypes.data_as(c_void_p),
       h.ctypes.data_as(c_void_p))
print "finished c"
for i in range(0,62):
    print(i)
    m2=map(lambda r: map( lambda e: doublify(e)[i],r), m)
    x = np.vectorize(lambda x: .25+int(x)/2.)(m2)
    plt.axis('on')
    imgplot = plt.imshow(x,vmin=0, vmax=1,interpolation='none',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
    plt.savefig("henonfigs/result_"+str(i)+".png")
    plt.axis('off')
    h2= map(lambda r: map(lambda e: doublify(e)[i],r),h)
    x2 = np.vectorize(lambda x: int(x)/1.)(h2)
    x2 = np.ma.masked_where(x2 <0.9,x2)
    imgplot3=plt.imshow(x3,vmin=0,vmax=1,interpolation='none',cmap=cm.Greys,extent= [xmin,xmax,ymin,ymax])
    imgplot2=plt.imshow(x2,vmin=0, vmax=1,interpolation='none',cmap=cm.Greens,extent=[xmin*mul,xmax*mul,ymin*mul,ymax*mul],alpha=0.6)
    plt.savefig("henondiffs/diff_"+str(i)+".png")
    plt.clf()
end=time.time()
print "total time: "+str(end-start)


