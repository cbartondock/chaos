import numpy as np
import os
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle
from ctypes import *
import time
from optparse import OptionParser

start=time.time()
parser = OptionParser()
parser.add_option('-g','--graph', action='store_true',dest='graph')
(options, args) = parser.parse_args()

#Required Libraries
libraries = open("../libraries.txt").readlines()
INV = CDLL(libraries[0].strip())
SAM = CDLL(libraries[2].strip())


#Choose a Map: (1=Standard,2=Pendulum Poincare Section,3=Henon,4=Saddle)
chmap = 4

#Topology: (1=Plane,2=Cylinder)
top = 2 if chmap == 2 else 1

#Iteration Properties
maxiter = 30
numper = 1

#Choose Grid size:
grid = 50

#Region Dimensions (make sure these are floats):
if chmap == 1:
    ymin = 0.
    ymax = 2.*np.pi     #Standard
    xmin = 0.
    xmax = 2.*np.pi
elif chmap == 2:
    ymin = -2.
    ymax = 4.        #Pendulum
    xmin = 0.
    xmax = 2.*np.pi

else:
    ymin=-2.
    ymax=2.      #Henon, Saddle
    xmin=-2.
    xmax=2.

deltax = (xmax-xmin)/grid
deltay = (ymax-ymin)/grid

#Property Summary
mapnames = ["Standard","Pendulum","Henon","Saddle"]
print("Running on {0} with box [({1},{2}),({3},{4})] with grid {5}, number per iterate {6}, and max iterates {7}".format(mapnames[chmap-1],xmin,ymin,xmax,ymax,grid,numper,maxiter))

params = np.savetxt("outputs/parameters.txt", np.array([int(grid), float(xmin),float(ymin),float(deltax),float(deltay)]))

m = np.ones((grid,grid),dtype="uint8")

print "started invariant4.c"
INV.calc_invariant(c_int32(maxiter),
       c_int32(maxiter),
       c_int32(numper),
       c_int32(grid),
       c_int32(grid),
       c_int32(chmap),
       c_int32(top),
	   c_longdouble(xmin),
	   c_longdouble(ymin),
	   c_longdouble(deltax),
       c_longdouble(deltay),
	   m.ctypes.data_as(c_void_p))
print "finished invariant4.c"


"""
#Plotting
w = np.vectorize(lambda x: x)(m)
plt.imshow(w,vmin=0, vmax=1,interpolation='none' if grid > 200 else 'nearest',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
plt.title("Invariant Set for the " + mapnames[chmap-1] + " Map")
if top == 1:
    plt.xlabel("x")
    plt.ylabel("y")
elif top==2:
    plt.xlabel("theta")
    plt.ylabel("omega")

plt.savefig("outputs/invariance_result.ps")
plt.clf()

#Saving
np.savetxt("outputs/imatrix.txt",w)
print("invariant saved")

#Graph Generation
class adj_element(Structure):
    pass
adj_element._fields_ = [("imageindex",c_int), ("domindex",c_int), ("imagenumber",c_int), ("next",POINTER(adj_element)), ("prev",POINTER(adj_element))]
class sparse_adjacency_matrix(Structure):
    _fields_=[("domnumber",c_int),("grid",c_int),("leastx",c_double),("leasty",c_double),("deltax",c_double),("deltay",c_double), ("adjacency_lists",POINTER(POINTER(adj_element)))]

sam_pointer = POINTER(sparse_adjacency_matrix)
adj_pointer = POINTER(adj_element)
initialize_matrix = SAM.initialize_sparse_matrix
initialize_matrix.restype = sam_pointer
free_matrix = SAM.free_matrix

if options.graph:
    sp = sam_pointer()
    sp = initialize_matrix(c_int(grid),c_int(numper),c_int(chmap),c_int(top),c_longdouble(xmin),c_longdouble(ymin),c_longdouble(deltax),c_longdouble(deltay),m.ctypes.data_as(c_void_p))
    nodes = {}
    current = adj_pointer()
    for i in range(0, sp.contents.domnumber):
        current = sp.contents.adjacency_lists[i].contents.next
        image=[]
        while current.contents.imageindex!=-1:
            image.append(current.contents.imageindex)
            current = current.contents.next
        nodes[current.contents.domindex] = image
    outfile = open("outputs/graph_save.p","wb")
    pickle.dump(nodes, outfile)
    outfile.close()
    free_matrix(sp)
    print "graph saved"

"""
