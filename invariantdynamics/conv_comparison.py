import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
import time

start = time.time()

hf = open("outputs/conv_hist.txt", "r")
histdata = hf.read().split(";")[:-1]
hf.close()

for hd in histdata:
    hd = hd.split('.')
    N = int(hd[0])
    print N
    numzeros = [int(i) for i in hd[1].split(', ')[:-1]]
    plt.hist(numzeros, np.arange(min(numzeros),max(numzeros),1)-1.5,alpha=0.5,label = "N={0}".format(N))
plt.legend(loc='upper right')
plt.xlabel('#zeros')
plt.ylabel('Number')
plt.savefig("outputs/multiple_conv_histograms.pdf")

""""
#Plotting Convergence Rate Histogram
entries = m.flatten()
plt.hist(entries,np.arange(min(entries),max(entries),1)-1.5)
plt.title("Convergence Rates Histogram (N={0})".format(totaltime))
plt.xlabel("#zeros")
plt.ylabel("N")
plt.savefig("outputs/convergence_histogram.ps")
plt.clf()
#Saving Convergence Matrix and histogram
with open("outputs/conv_matrix.p", "wb") as outfile:
    pickle.dump(m,outfile)
with open("outputs/conv_hist.txt", "a") as myfile:
        myfile.write(str((totaltime,list(entries))))
"""
end=time.time()
print(end-start)
