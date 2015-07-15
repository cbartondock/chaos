from matplotlib import pyplot as plt
import numpy as np
import math

f = open('histo.txt','r')
numx = []
numy = []
count = 0
for line in f:
    line = line.split()[:2]
    numx.append(int(line[0]))
    numy.append(int(line[1]))
    count += 1
numx.append(count)

plt.bar(numx[:-1],numy,width=1)
plt.xlim(min(numx), max(numx))
plt.show()

f.close()

