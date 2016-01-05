#Program for Distinguishing Quasiperiodic sets from Chaotic sets using WB averages

import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
from optparse import OptionParser
import time

start = time.time()
parser = OptionParser()
parser.add_option('-o','--opt', action='store_true',dest='optimized')
(options, args) = parser.parse_args()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
if options.optimized:
    QUAS=CDLL(liblist[14].strip())
else:
    QUAS=CDLL(liblist[4].strip())
libfile.close()

#Parameters

#Choose a Grid Size
grid=50

#Number of Functions Converging
fnum = 1

#Region Parameters
xmin=.7*np.pi
ymin=-.3*np.pi
xmax=1.3*np.pi
ymax=.3*np.pi
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

#Iteration Number
totaltime = 10000  # Map iterates

#Save Parameters
np.savetxt("outputs/qparameters.txt", np.array([int(grid), float(xmin),float(ymin),float(deltax),float(deltay)]))

m = -1*np.ones((grid,grid),dtype="int8")


print "started quasiperiodicity.c"
QUAS.convergence(c_int32(grid),
                 c_int32(grid),
                 c_int32(totaltime),
                 c_longdouble(xmin),
                 c_longdouble(ymin),
                 c_longdouble(deltax),
                 c_longdouble(deltay),
		         c_int32(fnum),
                 m.ctypes.data_as(c_void_p))
print "finished quasiperiodicity.c"

#Plotting Convergence Spatially
w = np.vectorize(lambda x: x)(m)
print(m)
plt.imshow(w,vmin=0,interpolation='nearest',cmap=cm.jet,extent=[xmin,xmax,ymin,ymax],origin='lower')
cbar = plt.colorbar()
#plt.title("Rates of Birkhoff Convergence in the Standard Map (N={0})".format(totaltime))
cbar.set_label("#zeros")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("outputs/convergence_result_t{0}_g{1}_xs{2}_ys{3}_xb{4}_yb{5}.pdf".format(totaltime,grid,round(xmin,2),round(ymin,2),round(xmax,2),round(ymax,2)))
plt.clf()

#Plotting Convergence Rate Histogram
entries = m.flatten()
plt.hist(entries,np.arange(min(entries),max(entries),1)-1.5, alpha=0.5)
#plt.title("Convergence Rates Histogram (N={0})".format(totaltime))
plt.xlabel("#zeros")
plt.ylabel("N")
plt.savefig("outputs/convergence_histogram_t{0}_g{1}_xs{2}_ys{3}_xb{4}_yb{5}.pdf".format(totaltime,grid,round(xmin,2),round(ymin,2),round(xmax,2),round(ymax,2)))
plt.clf()
#Saving Convergence Matrix and histogram
with open("outputs/conv_matrix.p", "wb") as outfile:
    pickle.dump(m,outfile)
with open("outputs/conv_hist.txt", "a") as myfile:
    myfile.write(str(totaltime)+"."+str(list(entries))[1:-2]+";")

end=time.time()
print(end-start)
