import time
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




xmin=-2
xmax=2
ymin=-.88
ymax=.23
eps=.1

def alg(xmin,xmax,ymin,ymax,eps, n,numpoints,smap):
    if n==0:
        return [[xmin,xmax,ymin,ymax]]
    print(n)
    pointsup = pickn(numpoints,xmin,xmax,ymax,ymax+eps)
    pointsdown = pickn(numpoints,xmin,xmax,ymin-eps,ymin)
    pointsright = pickn(numpoints,xmax,xmax+eps,ymin,ymax)
    pointsleft = pickn(numpoints, xmin-eps,xmin,ymin,ymax)
    
    pointsup = [smap(p) for p in pointsup]
    pointsdown = [smap(p) for p in pointsdown]
    pointsright = [smap(p) for p in pointsright]
    pointsleft = [smap(p) for p in pointsleft]
    
    vdegree=0

    for p in pointsup:
        if p[1]<ymax and p[1]>=ymin and p[0]<xmax and p[0]>=xmin:
            vdegree+=1
            break
    for p in pointsdown:
        if p[1]<ymax and p[1]>=ymin and p[0]<xmax and p[0]>=xmin:
            vdegree+=1
            break
    for p in pointsright:
        if p[1]<ymax and p[1]>=ymin and p[0]<xmax and p[0]>=xmin:
            vdegree+=1
            break
    for p in pointsleft:
        if p[1]<ymax and p[1]>=ymin and p[0]<xmax and p[0]>=xmin:
            vdegree+=1
            break

    if vdegree >=2:
        return alg(xmin,0.5*(xmin+xmax),ymin,0.5*(ymin+ymax),eps,n-1,numpoints,smap) + alg(.5*(xmin+xmax),xmax,ymin, 0.5*(ymin+ymax),eps,n-1,numpoints,smap)+alg(xmin,0.5*(xmax+xmax),0.5*(ymin+ymax),ymax,eps,n-1,numpoints,smap) + alg(0.5*(xmin+xmax),xmax,0.5*(ymax+ymin),ymax,eps,n-1,numpoints,smap)
    else:
        return []



result=alg(xmin,xmax,ymin,ymax,.1,5,10,testmap)
print(result)

plt.axes()
domain = plt.Rectangle((xmin,ymin),xmax-xmin,ymax-ymin,fc='r')
plt.gca().add_patch(domain)
for r in result:
    rectangle = plt.Rectangle((r[0],r[2]),r[1]-r[0],r[3]-r[2],fc='b')
    plt.gca().add_patch(rectangle)
circle = plt.Circle((0,0),.01,fc='y')
plt.gca().add_patch(circle)
plt.axis('scaled')
plt.savefig('alg1_result.png')




