# This file was *autogenerated* from the file smap4.sage
from sage.all_cmdline import *   # import sage library
_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2p5 = RealNumber('2.5'); _sage_const_5 = Integer(5); _sage_const_100 = Integer(100); _sage_const_1p4 = RealNumber('1.4'); _sage_const_p5 = RealNumber('.5'); _sage_const_p2 = RealNumber('.2'); _sage_const_p3 = RealNumber('.3')
"""Dash's Algorithm"""

import time
from collections import Counter
from random import random as rand
import matplotlib.pyplot as plt
from sets import Set

#modulus
def rmod(x, y):
    result= x - int(x/y)*y
    return result if result >= _sage_const_0  else result + y

#henon map generator
def hmap_g(a,b):
    def hmap(p):
        x= N(a-p[_sage_const_0 ]**_sage_const_2 +b*p[_sage_const_1 ])
        y= N(p[_sage_const_0 ])
        return [x,y]
    return hmap


#standard map generator
def smap_g(rho):
    def smap(p):
        x = N(rmod(p[_sage_const_0 ]+p[_sage_const_1 ],_sage_const_2 *pi))
        y = N(rmod(p[_sage_const_1 ] + rho*sin(p[_sage_const_0 ]+p[_sage_const_1 ]),_sage_const_2 *pi))
        return [x,y]
    return smap

#point picker
def pickpoint(xmin,xmax,ymin,ymax):
    return [rand()*(xmax-xmin)+xmin,rand()*(ymax-ymin)+ymin]
def pickn(n,xmin,xmax,ymin,ymax):
    result=[]
    for i in range(_sage_const_0 ,n):
        result.append(pickpoint(xmin,xmax,ymin,ymax))
    return result

#maps
myhmap = hmap_g(_sage_const_1p4 ,_sage_const_p3 )
mysmap = smap_g(_sage_const_1 )
testmap = lambda p: [_sage_const_p5 *p[_sage_const_0 ],_sage_const_3 *p[_sage_const_1 ]]
testmap2 = lambda p: [_sage_const_5 *p[_sage_const_0 ],_sage_const_p3 *p[_sage_const_1 ]]
rotmap = lambda theta: lambda p: [cos(theta)*p[_sage_const_0 ]-sin(theta)*p[_sage_const_1 ],sin(theta)*p[_sage_const_0 ]+cos(theta)*p[_sage_const_1 ]]
"""xmin=0
xmax=2*N(pi)
ymin=0
ymax=2*N(pi)
"""
xmin=-_sage_const_2p5 
xmax=_sage_const_2p5 
ymin=-_sage_const_2p5 
ymax=_sage_const_2p5 



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
                if p[_sage_const_0 ] < xmax and p[_sage_const_0 ]>=xmin and p[_sage_const_1 ] < ymax and p[_sage_const_1 ]>=ymin:
                    boxes.add( (floor(p[_sage_const_0 ]/deltax)*deltax,floor(p[_sage_const_1 ]/deltay)*deltay) )
            x+=deltax
        y+=deltay
    s=False
    while n>_sage_const_0 :
        print(n)
        n-=_sage_const_1 
        newboxes=Set([])
        for box in boxes:
            points = pickn(numpoints,box[_sage_const_0 ],box[_sage_const_0 ]+deltax, box[_sage_const_1 ],box[_sage_const_1 ]+deltay)
            points = [smap(p) for p in points]
            for p in points:
                print p
                if p[_sage_const_0 ] < xmax and p[_sage_const_0 ]>=xmin and p[_sage_const_1 ] < ymax and p[_sage_const_1 ]>=ymin:
                    newboxes.add( (floor(p[_sage_const_0 ]/deltax)*deltax,floor(p[_sage_const_1 ]/deltay)*deltay) )
        if newboxes==boxes:
            s=True
            break
        boxes=newboxes
    for box in boxes:
        rectangle=plt.Rectangle((box[_sage_const_0 ],box[_sage_const_1 ]),deltax,deltay,fc='b',ec='none')
        plt.gca().add_patch(rectangle)
    print "successful" if s else "failure"
    plt.axis('scaled')
    plt.ylim([ymin,ymax])
    plt.xlim([xmin,xmax])
    plt.savefig("alg4_result.png")
    


    
alg4(xmin,xmax,ymin,ymax,_sage_const_100 ,rotmap(_sage_const_p2 ),_sage_const_100 ,_sage_const_5 )
    

