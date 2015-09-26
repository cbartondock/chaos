
#Program for Distinguishing Quasiperiodic sets from Chaotic sets using WB averages

import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import cPickle as pickle;
from ctypes import *
import time

start = time.time()

with f as open("text_quasi_conv_t"):
    lines = f.read().split("\n")



#Plotting Convergence Spatially
w = np.vectorize(lambda x: x)(m)
print(m)
plt.imshow(w,vmin=0,interpolation='nearest',cmap=cm.Blues,extent=[xmin,xmax,6.28-ymax,6.28-ymin])
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
