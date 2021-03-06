from optparse import OptionParser
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import time

parser = OptionParser()
parser.add_option('-f','--file',dest='filename')
(options, args) = parser.parse_args()

grid = int(options.filename.split("_")[4][1:])
time = int(options.filename.split("_")[3][1:])
xmin = float(options.filename.split("_")[5][2:])
ymin = float(options.filename.split("_")[6][2:])
xmax = float(options.filename.split("_")[7][2:])
ymax = float(options.filename.split("_")[8].split(".txt")[0][2:])
name = "outputs/text_quasi_conv_t{0}_g{1}_xs0.00_ys0.00_xb6.28_yb6.28.txt".format(time,grid)
f = open(options.filename,"r")
data= []
for line in f.readlines():
    data.append([float(item.split(": ")[1]) for item in line.split(', ')[2:]])
m = []
for i in range(0, grid):
    m.append([])
    for j in range(0, grid):
        m[i].append(data[grid*i+j][2])

plt.imshow(m,vmin=0,interpolation='nearest',cmap=cm.jet,extent=[xmin,xmax,ymin,ymax],origin = 'lower')
cbar = plt.colorbar()
cbar.set_label("#zeros")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("outputs/pconvergence_result_t{0}_g{1}_xs{2}_ys{3}_xb{4}_yb{5}.pdf".format(time,grid,round(xmin,2),round(ymin,2),round(xmax,2),round(ymax,2)))
plt.clf()
f.close()

""""
start = time.time()

libfile = open("../libraries.txt")
liblist = libfile.readlines()
QUAS=CDLL(liblist[4].strip())
libfile.close()

#Parameters

#Choose a Grid Size
grid=2000

#Number of Functions Converging
fnum = 1

#Region Parameters
xmin=0*np.pi#0*np.pi
ymin=0*np.pi#-1*np.pi
xmax=2*np.pi#2*np.pi
ymax=2*np.pi#1*np.pi
deltax=(xmax-xmin)/grid
deltay=(ymax-ymin)/grid

#Iteration Number
totaltime = 200  # Map iterates

#Save Parameters
np.savetxt("outputs/qparameters.txt", np.array([int(grid), float(xmin),float(ymin),float(deltax),float(deltay)]))

m = np.ones((grid,grid),dtype="uint8")


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
plt.imshow(w,vmin=0,interpolation='nearest',cmap=cm.jet,extent=[xmin,xmax,6.28-ymax,6.28-ymin])
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
"""
