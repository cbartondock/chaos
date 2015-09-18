import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from ctypes import *
import time

start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
WCONV = CDLL(liblist[12].strip())
libfile.close()

#Trajectory Origin
x = 3.14
y = 0.1

totaltime = 1000 # Map iterates per window
num = 1000 #number of windows
cd = np.zeros(num)

print "started lyapunovconv.c"
WCONV.convergence(
        c_double(x),
        c_double(y),
        c_int32(totaltime),
        c_int32(num),
        c_int32(1),
        cd.ctypes.data_as(c_void_p))
print "finished lyapunovconv.c"
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.set_yscale('log')
ax.plot(cd[:-2])

plt.savefig("outputs/rates/wq_birkrate.png")
plt.clf()
end = time.time()
print(end-start)
