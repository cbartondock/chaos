import numpy
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle
from ctypes import *
import time
from optparse import OptionParser

def rmod(x, y):
    result = x - int(x/y)*y
    return result if result >= 0 else result + y



start=time.time()
parser = OptionParser()
parser.add_option('-g','--graph', action='store_true',dest='graph')
(options, args) = parser.parse_args()

INV = CDLL('/Users/chrisdock/Documents/chaos/ccodefromsage/invariantdynamics/invariant4.dylib')
RK4 = CDLL('/Users/chrisdock/Documents/chaos/ccodefromsage/rk4/rk4.dylib')
SAM = CDLL('/Users/chrisdock/Documents/chaos/ccodefromsage/sparse_matrix_table/smtable.dylib')

maxiter = 30
grid = 1000
print "grid: "+str(grid)+", maxiter: "+str(maxiter)
"""
ymin=-2.5
ymax=2.5
xmin=-2.5
xmax=2.5
"""

ymin = -2
ymax = 4
xmin = 0
xmax = 2*numpy.pi

"""
ymin = 0
ymax = 2*numpy.pi
xmin = 0
xmax = 2*numpy.pi
"""

deltax = (xmax-xmin)/grid
deltay = (ymax-ymin)/grid

params = numpy.savetxt("outputs/parameters.txt", numpy.array([int(grid), float(xmin),float(ymin),float(deltax),float(deltay)]))



m = matrix.ones(grid,grid)
m = m.numpy('uint8')
print "started invariant.c"
INV.calc_invariant(c_int32(maxiter),
       c_int32(maxiter),
       c_int32(grid),
       c_int32(grid),
	   c_double(xmin),
	   c_double(ymin),
	   c_double(deltax),
       c_double(deltay),
	   m.ctypes.data_as(c_void_p))
print "finished invariant.c"
w = numpy.vectorize(lambda x: x)(m)
imgplot = plt.imshow(w,vmin=0, vmax=1,interpolation='none' if grid > 200 else 'nearest',cmap=cm.Blues,extent=[xmin,xmax,ymin,ymax])
plt.savefig("outputs/invariance_result.ps")
numpy.savetxt("outputs/imatrix.txt",w)
print("invariant saved")

plt.clf()

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
    sp = initialize_matrix(c_int(grid),c_double(xmin),c_double(ymin),c_double(deltax),c_double(deltay),m.ctypes.data_as(c_void_p))
    nodes = {}
    current = adj_pointer()
    for i in range(0, sp.contents.domnumber):
        current = sp.contents.adjacency_lists[i].contents.next
        image=[]
        while current.contents.imageindex!=-1:
            image.append(current.contents.imageindex)
            current = current.contents.next
        nodes[current.contents.domindex] = image
    outfile = open("outputs/save.p","wb")
    pickle.dump(nodes, outfile)
    outfile.close()
    free_matrix(sp)
    print "graph saved"


    

