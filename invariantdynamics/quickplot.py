import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
from optparse import OptionParser
import time



parser = OptionParser()
parser.add_option('-f','--file',dest='filename')
(options, args) = parser.parse_args()

with open(options.filename) as f:
    lines = f.readlines();
    point_list= [(float(line.split(":")[3].split(',')[0].strip()),float(line.split(":")[4].split(",")[0].strip()), line) for line in lines]
    xs = [p[0] for p in point_list]
    ys =  [-1*p[1] for p in point_list]
    plt.scatter(xs,ys,s=10)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
end = time.time()
