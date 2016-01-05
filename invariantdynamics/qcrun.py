

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
ax = 0#np.pi
bx = 2*np.pi#.84*1.55*np.pi
ay = 0#0
by = 0#.84

#Iteration Number
totaltime = 1000000   # Map iterates
numpoints = 1000 #sampling of curve
totalshifts = 1
alldata = [[0 for i in range(totalshifts)] for j in range(numpoints)]
#Save Parameters
for preshift in range(totalshifts):
    m = np.zeros((numpoints,2),dtype="float64")
    print preshift
    #print "started quasiperiodicity.c"
    CQUAS.curve_convergence(c_longdouble(ax),
                 c_longdouble(ay),
                 c_longdouble(bx),
                 c_longdouble(by),
                 c_int32(numpoints),
                 c_int32(totaltime),
		         c_int32(fnum),
                 c_int32(0),
                 m.ctypes.data_as(c_void_p))
    #print "finished quasiperiodicity.c"
    m = np.transpose(m)

    plt.plot(m[0],m[1],'-', color='b')

    for i in range(len(m[1])):
        alldata[i][preshift]=m[1][i]
averages = [sum(data)/float(totalshifts) for data in alldata]
#print averages
#variances = [(sum([(alldata[i][j] - averages[i])**2 for j in range(totalshifts)])**.5)/float(totalshifts) for i in range(numpoints)]
#print max(variances)
print "finished all"


plt.xlabel('parameter')
plt.ylabel('#zeros')
plt.savefig('outputs/quasi_curve_result_t{0}_np{1}_ax{2}_ay{3}_bx{4}_by{5}.pdf'.format(totaltime,numpoints,round(ax,2),round(ay,2),round(bx,2),round(by,2)))
plt.clf()
#plt.plot(m[0],variances,'-')
#plt.xlabel('parameter')
#plt.ylabel('Average Deviation')
#plt.savefig('outputs/quasi_variance_result.pdf')
end = time.time()
print(end-start)

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
