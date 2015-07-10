"""Dash's Algorithm"""

import time
from collections import Counter
from random import random as rand
import matplotlib.pyplot as plt
from sets import Set

#modulus
def rmod(x, y):
    result= x - int(x/y)*y
    return result if result >= 0 else result + y

#henon map generator
def hmap_g(a,b):
    def hmap(p):
        x= N(a-p[0]^2+b*p[1])
        y= N(p[0])
        return [x,y]
    return hmap


#standard map generator
def smap_g(rho):
    def smap(p):
        x = N(rmod(p[0]+p[1],2*pi))
        y = N(rmod(p[1] + rho*sin(p[0]+p[1]),2*pi))
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
testmap2 = lambda p: [5*p[0],.3*p[1]]
rotmap = lambda theta: lambda p: [cos(theta)*p[0]-sin(theta)*p[1],sin(theta)*p[0]+cos(theta)*p[1]]
"""xmin=0
xmax=2*N(pi)
ymin=0
ymax=2*N(pi)
"""
xmin=-2.5
xmax=2.5
ymin=-2.5
ymax=2.5



def alg4(xmin,xmax,ymin,ymax,grid,smap,numpoints,n):
    plt.axes()
    deltax= (xmax-xmin)/grid
    deltay= (ymax-ymin)/grid
    
    boxes = Set([])
    y=ymin
    while y < ymax:
        x=xmin
        while x < xmax:
            print(y)
            rectangle= plt.Rectangle((x,y),deltax,deltay,fc='g',ec='none')
            plt.gca().add_patch(rectangle)
            points = pickn(numpoints,x,x+deltax,y,y+deltay)
            points = [smap(p) for p in points]
            for p in points:
                if p[0] < xmax and p[0]>=xmin and p[1] < ymax and p[1]>=ymin:
                    boxes.add( (floor(p[0]/deltax)*deltax,floor(p[1]/deltay)*deltay) )
            x+=deltax
        y+=deltay
    s=False
    while n>0:
        print(n)
        n-=1
        newboxes=Set([])
        for box in boxes:
            points = pickn(numpoints,box[0],box[0]+deltax, box[1],box[1]+deltay)
            points = [smap(p) for p in points]
            for p in points:
                print p
                if p[0] < xmax and p[0]>=xmin and p[1] < ymax and p[1]>=ymin:
                    newboxes.add( (floor(p[0]/deltax)*deltax,floor(p[1]/deltay)*deltay) )
        if newboxes==boxes:
            s=True
            break
        boxes=newboxes
    for box in boxes:
        rectangle=plt.Rectangle((box[0],box[1]),deltax,deltay,fc='b',ec='none')
        plt.gca().add_patch(rectangle)
    print "successful" if s else "failure"
    plt.axis('scaled')
    plt.ylim([ymin,ymax])
    plt.xlim([xmin,xmax])
    plt.savefig("alg4_result.png")
    


    
alg4(xmin,xmax,ymin,ymax,100,rotmap(.2),100,5)
    

