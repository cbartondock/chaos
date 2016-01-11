""" This is the governor program for birk_curve.c. Its inputs are two points, the number of sample points
on the line between those two points, and the amount of iteration time used for each initial condition. Its output is a birkhoff average associated to each point on the parametrized segment.
"""

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

#Plotting
fig = plt.figure()
fig.suptitle(r'Birkhoff Average of $\phi$ along parametrized segment'.format(m),fontsize=14)
axes = fig.add_subplot(111)
fig.subplots_adjust(top=.85)
axes.plot(m[0],m[1],'-')
axes.set_xlabel(r'$\lambda: [0,1]\mapsto (1-t)\cdot({0},{1}) + t\cdot({2},{3})$'.format(ax,ay,bx,by))
axes.set_ylabel(r'$\langle\phi(\lambda)\rangle$')
plt.savefig('outputs/birk_curve_result.pdf')
plt.clf()



"""
    fig = plt.figure()
    fig.suptitle(r'Analysis of Rational Rotation Numbers $\frac{{p}}{{q}}$ with $p+q<{0}$'.format(m),fontsize=14,fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=.85)

    ax.plot(lambdalistq,nzlist,'-', color='b',linewidth=.8)
    ax.plot(lambdalistr,rholist,'-', color = 'r',linewidth=.8)

    pairlist = zip(rholist,lambdalistr)
    k=0
    for r in rationals:
        rat = float(r[0])/float(r[1])
        mindiff = 10
        minpair = (0,0)
        for pair in pairlist:
            if abs(pair[0]-rat) < mindiff:
                mindiff = abs(pair[0]-rat)
                minpair = pair
        if abs(minpair[0]-rat) < .005:
            ax.plot([minpair[1],minpair[1]],[minpair[0]-.005,minpair[0]+.01],'-', color = 'k', linewidth =.4)
            ax.text(minpair[1]-.0035,minpair[0]+.013+.011*k,r'$\frac{{{0}}}{{{1}}}$'.format(r[0],r[1]),fontsize=5)
        k=1-k

    ax.set_xlabel(r'$\lambda: [0,1]\mapsto (1-t)\cdot({0},{1}) + t\cdot({2},{3})$'.format(axq,ayq,bxq,byq))
    ax.set_ylabel(r'$\#zeros$ (blue), $\rho$ (red)')

    out = qfilename.split("/")[0]+"/"+"pquasirot_curve_result_t"+qfilename.split("/")[1].split("_t")[1]
    out = out[:-4]+"_r{0}".format(m)+".pdf"
    plt.savefig(out)
    plt.clf()
    fq.close()
    fr.close()

"""
