"""Dash's Algorithm"""

import time
from collections import Counter
from random import random as rand
import matplotlib.pyplot as plt
from sets import Set

#modulus
def rmod(x, y):
    result= x - int(x/y)*y;
    return result if result > 0 else y + result    

#standard map generator
def smap_g(rho):
    def smap(p):
        x = rmod(N(p[0]+p[1]),2*N(pi))
        y = rmod(N(p[1] + rho*sin(p[0]+p[1])),2*N(pi))
        return [x,y]
    return smap


mysmap = smap_g(1)
trajx=[]
trajy=[]
point=(0.1,0)
for i in range(0,10000):
    print(i)
    trajx.append(point[0])
    trajy.append(point[1])
    point=mysmap(point)
plt.xlim([0,6.28])
plt.ylim([0,6.28])
lp=plt.scatter(trajx,trajy)
plt.savefig("trajectory.png")
