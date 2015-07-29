import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
import cPickle as pickle
import time
from math import floor

start = time.time()

params = list(np.loadtxt("outputs/qparameters.txt"))
print(params)
infile = open("outputs/conv_matrix.p")
quasimatrix = pickle.load(infile)
infile.close()

pointsfile = open("outputs/stickypoints.txt")
stickymatrix = np.zeros((params[0],params[0]))
k=1
for pointsstring in pointsfile.read().split(';'):
    k+=1
    pointsstring = pointsstring[:-2]
    points = pointsstring.split("(")[1:]
    points=[point[:-3].split(", ") for point in points]
    points = [[float(coord) for coord in point] for point in points]
    print(len(points))
    for w in range(int(.4*len(points)),int(.62*len(points))):
        point = points[w]
        y = point[0]
        x = point[1]
        i = floor((x-params[1])/params[3])
        j = floor((y-params[2])/params[4])
        stickymatrix[i][j]=k
        stickymatrix[i-1][j]=k
        stickymatrix[i+1][j]=k
        stickymatrix[i][j-1]=k
        stickymatrix[i][j+1]=k
pointsfile.close()
stickymatrix = np.ma.masked_where(stickymatrix < 1,stickymatrix)
plt.imshow(quasimatrix,interpolation='nearest',cmap=cm.Blues,extent=[params[1],params[1]+params[0]*params[3],params[2],params[2]+params[0]*params[4]])
cbar = plt.colorbar()
cbar.set_label("#zeros")
plt.xlabel("x")
plt.ylabel("y")
#plt.title("Sticky Points of a Chaotic Trajectory")
plt.imshow(stickymatrix,interpolation='sinc',cmap=cm.brg,extent=[params[1],params[1]+params[0]*params[3],params[2],params[2]+params[0]*params[4]])
red_patch = mpatches.Patch(color='#f90504', label='Sticky Region II')
green_patch = mpatches.Patch(color='#00fd39', label='Sticky Region III')
blue_patch = mpatches.Patch(color='#0025f6', label='Sticky region I')
fontP = FontProperties()
fontP.set_size('small')
plt.legend(handles=[blue_patch,red_patch,green_patch],loc = "upper right",bbox_to_anchor = (0.17, 1),prop=fontP)
plt.savefig("outputs/sticky_result.png")
plt.clf()


end=time.time()
print(end-start)
