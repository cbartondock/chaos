"""
This program plots the rotation number against the convergence rate of the birkhoff average.
It takes two datafiles as inputs, a convergence file and a rotation # file. It is run as follows:
    python rotquasi_plotter text_quasi_curve[].txt text_rot_[].txt
"""

from optparse import OptionParser
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from matplotlib import rc
import time
import simplest_rationals as sr

ms = [15]
rational_analysis=False

plt.rc('text', usetex=True)
for m in ms:

    rationals = sr.simplest_rationals_between(5/38.+.00001,1/5.,m)
    rationals.sort(key = lambda x: float(x[0])/float(x[1]))
    print(rationals)

    parser = OptionParser()
    ( args) = parser.parse_args()
    qfilename = str(args[1][0])
    rfilename = str(args[1][1])
    print qfilename
    print rfilename

    numpointsq = int(qfilename.split("_")[4][2:])
    timeq = int(qfilename.split("_")[3][1:])
    axq = float(qfilename.split("_")[5][2:])
    ayq = float(qfilename.split("_")[6][2:])
    bxq = float(qfilename.split("_")[7][2:])
    byq = float(qfilename.split("_")[8].split(".txt")[0][2:])

    numpointsr = int(rfilename.split("_")[4][2:])
    timer = int(rfilename.split("_")[3][1:])
    axr = float(rfilename.split("_")[5][2:])
    ayr = float(rfilename.split("_")[6][2:])
    bxr = float(rfilename.split("_")[7][2:])
    byr = float(rfilename.split("_")[8].split(".txt")[0][2:])

    if axr!=axq or ayr!=ayq or bxr!=bxq or byr!=byq:
        print "file mismatch"
        exit()

    fq = open(qfilename, "r")
    nzlist= []
    for line in fq.readlines():
        nzlist.append(float(line.split("numzeros:")[1].strip()))
    lambdalistq = [float(i)/float(numpointsq) for i in range(0,numpointsq)]
    maxnz = max([nz for nz in nzlist if nz <1000])

    fr = open(rfilename,"r")
    rholist = []
    for line in fr.readlines():
        rholist.append(float(line.split("rho:")[1].strip()))
    lambdalistr= [float(i)/float(numpointsr) for i in range(0,numpointsr)]
    maxrho = max([rho for rho in rholist if rho<1000])
    print(maxnz)
    nzlist = [maxrho*nz/maxnz for nz in nzlist]

    fig = plt.figure()
    fig.suptitle(r'Analysis of Rational Rotation Numbers $\frac{{p}}{{q}}$ with $p+q<{0}$'.format(m),fontsize=14,fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=.85)

    ax.plot(lambdalistq,nzlist,'-', color='b',linewidth=.8)
    ax.plot(lambdalistr,rholist,'-', color = 'r',linewidth=.8)
    if rational_analysis:
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
