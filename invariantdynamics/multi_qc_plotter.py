from optparse import OptionParser
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import time
import matplotlib.patches as mpatches

parser = OptionParser()
(args) = parser.parse_args()
f1_name = str(args[1][0])
f2_name = str(args[1][1])
f1 = open(f1_name)
f2 = open(f2_name)


name = "outputs/text_quasi_curve_t20000_np1000_ax0.00_ay0.00_bx6.28_by0.00.txt"

numpoints1 = int(f1_name.split("_")[4][2:])
time1 = int(f1_name.split("_")[3][1:])
ax1 = float(f1_name.split("_")[5][2:])
ay1 = float(f1_name.split("_")[6][2:])
bx1 = float(f1_name.split("_")[7][2:])
by1 = float(f1_name.split("_")[8].split(".txt")[0][2:])

numpoints2 = int(f2_name.split("_")[4][2:])
time2 = int(f2_name.split("_")[3][1:])
ax2 = float(f2_name.split("_")[5][2:])
ay2 = float(f2_name.split("_")[6][2:])
bx2 = float(f2_name.split("_")[7][2:])
by2 = float(f2_name.split("_")[8].split(".txt")[0][2:])


nzlist1= []
nzlist2= []
for line in f1.readlines():
    nzlist1.append(float(line.split("numzeros:")[1].strip()))
for line in f2.readlines():
    nzlist2.append(float(line.split("numzeros:")[1].strip()))

lambdalist = [float(i)/float(numpoints1) for i in range(0,numpoints1)]
plt.plot(lambdalist,nzlist1,'-', color='g')
plt.plot(lambdalist,nzlist2,'-', color='r')
gpatch = mpatches.Patch(color='green',label=r'$N=10^4$')
rpatch = mpatches.Patch(color='red',label=r'$N=10^6$')
plt.legend(handles=[gpatch,rpatch])
plt.xlabel("parameter")
plt.ylabel("numzeros")
out ="pquasi_curve_combined.pdf"
plt.savefig(out)
plt.clf()
f1.close()
f2.close()

