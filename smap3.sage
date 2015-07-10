import time
from collections import Counter
from random import random as rand
import matplotlib.pyplot as plt



#modulus
def rmod(x, y):
    return x - int(x/y)*y

#henon map generator
def hmap_g(a,b):
    def hmap(p):
        x= a-p[0]^2+b*p[1]
        y= p[0]
        return [x,y]
    return hmap


#standard map generator
def smap_g(rho):
    def smap(p):
        x = rmod(p[0]+p[1],2*N(pi))
        y = p[1] + rho*sin(p[0]+p[1])
        return [x,y]
    return smap

#point picker
def pickpoint(xmin,xmax,ymin,ymax):
    return [rand()*(xmax-xmin)+xmin,rand()*(ymax-ymin)+ymin]
def pickn(n,xmin,xmax,ymin,ymax):
    result=[]
    for i in range(0,n):
        result.append(pickpoint(xmin,xmax,ymin,ymax))
    return result

#maps
myhmap = hmap_g(1.4,.3)
mysmap = smap_g(1)
testmap = lambda p: [.5*p[0],3*p[1]]
testmap2 = lambda p: [p[0]+sin(3.14*(p[0]+1)),.3*p[1]]


xmin=-2.5
xmax=2.5
ymin=-2.5
ymax=2.5


def alg3(xmin,xmax,ymin,ymax,grid,smap,n,numpoints):
    plt.axes()
    deltax= (xmax-xmin)/grid
    deltay= (ymax-ymin)/grid
    points = pickn(numpoints, xmin, xmax, ymin, ymax)
    for i in range(0,n):
        points = [smap(p) for p in points]
    y=ymin
    recs=[]
    i=0
    cnt = Counter()
    for p in points:
        cnt[(floor(p[0]/deltax),floor(p[1]/deltay))]+=1
    while y<ymax:
        x=xmin
        while x<xmax:
            print(i)
            i+=1
            rectangle= plt.Rectangle((x,y),deltax,deltay,fc='g')
            plt.gca().add_patch(rectangle)
            if cnt[(floor(x/deltax),floor(y/deltay))]>0:
                rectangle= plt.Rectangle((x,y),deltax,deltay,fc='b')
                plt.gca().add_patch(rectangle)
            x+=deltax
        y+=deltay
    plt.axis('scaled')
    plt.ylim([-2,2])
    plt.xlim([-2,2])
    plt.savefig("alg3_result.png")
print "calling"
alg3(xmin,xmax,ymin,ymax,50.,mysmap,8,100000)
