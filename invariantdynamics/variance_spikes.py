


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
    rho = float(line.split("rho:")[1].strip())
    if rho==rho:
        rholist.append(float(line.split("rho:")[1].strip()))
print rholist
lambdalist = [float(i)/float(numpoints) for i in range(numpoints-len(rholist),numpoints)]
deg=1
def exclude(l):
    print(len(l))
    return [l[i] for i in range(0,len(l)) if (i<250 or i>350)]
print(exclude(lambdalist))
cfs = [-0.00073926, 0.1934058]
#cfs = np.polyfit(exclude(lambdalist), exclude(rholist), deg)
print cfs
fitlist = [sum([cfs[i]*(x**(deg-i)) for i in range(0,deg+1)]) for x in lambdalist]
vardict = [(fitlist[i], abs(fitlist[i]-rholist[i])) for i in range(0,len(rholist))]
print max(vardict, key = lambda p: p[1])
plt.plot(lambdalist, [(fitlist[i]-rholist[i])**2 for i in range(0,len(rholist))], '-',color='g')
plt.xlabel("parameter")
plt.ylabel("variance")
out = options.filename.split("/")[0]+"/"+"prot_var_t"+options.filename.split("/")[1].split("_t")[1]
out = out[:-3]+"pdf"
plt.savefig(out)
plt.clf()
plt.plot(lambdalist,fitlist,'-',color='g')
plt.plot(lambdalist,rholist,'-',color='r')
plt.show()
plt.clf()
f.close()

