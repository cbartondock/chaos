from optparse import OptionParser
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import time

parser = OptionParser()
parser.add_option('-f','--file',dest='filename')
(options, args) = parser.parse_args()
name = "outputs/text_rotation_curve_t1000_np100_ax3.14_ay0.00_bx4.08_by0.00.txt"

numpoints = int(options.filename.split("_")[4][2:])
time = int(options.filename.split("_")[3][1:])
ax = float(options.filename.split("_")[5][2:])
ay = float(options.filename.split("_")[6][2:])
bx = float(options.filename.split("_")[7][2:])
by = float(options.filename.split("_")[8].split(".txt")[0][2:])

f = open(options.filename, "r")
rholist= []
for line in f.readlines():
    rholist.append(float(line.split("rho:")[1].strip()))
lambdalist = [float(i)/float(numpoints) for i in range(0,numpoints)]
plt.plot(lambdalist,rholist,'-', color='b')
plt.xlabel("parameter")
plt.ylabel("rotation num")
out = options.filename.split("/")[0]+"/"+"prot_curve_result_t"+options.filename.split("/")[1].split("_t")[1]
out = out[:-3]+"pdf"
plt.savefig(out)
plt.clf()
f.close()

