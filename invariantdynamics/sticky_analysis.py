import numpy as np
import matplotlib.pylab as plt
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


pointsstring = pointsfile.read()[:-2]
points = pointsstring.split("(")[1:]
points = [point[:-3].split(", ") for point in points]
points = [[float(coord) for coord in point] for point in points]
pointsfile.close()
print(len(points))
for w in range(int(.4*len(points)),int(.6*len(points))):
    point = points[w]
    y = point[0]
    x = point[1]
    i = floor((x-params[1])/params[3])
    j = floor((y-params[2])/params[4])
    stickymatrix[i][j]=1

stickymatrix = np.ma.masked_where(stickymatrix < 1,stickymatrix)
plt.imshow(quasimatrix,interpolation='nearest',cmap=cm.Blues,extent=[params[1],params[1]+params[0]*params[3],params[2],params[2]+params[0]*params[4]])
plt.colorbar()
plt.imshow(stickymatrix,interpolation='nearest',cmap=cm.hsv,extent=[params[1],params[1]+params[0]*params[3],params[2],params[2]+params[0]*params[4]])

plt.savefig("outputs/sticky_result.ps")
plt.clf()


end=time.time()
print(end-start)
