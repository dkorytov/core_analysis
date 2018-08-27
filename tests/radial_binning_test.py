#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt
from random import gauss
import sys


def radial_weights(radius, r_bins):
    weights = np.zeros(len(r_bins))
    if(radius<r_bins[0]):
        return weights;
    if(radius != 0):
        start_angle = np.arcsin(r_bins[0]/radius)
    else:
        start_angle = 0
    #end_angle = pi/2
    for i in range(0,len(radial_bins)-1):
        if(radius<r_bins[i+1]):
            weights[i] = np.cos(start_angle)
            break
        else:
            end_angle = np.arcsin(r_bins[i+1]/radius)
            area = -np.cos(end_angle) + np.cos(start_angle)
            weights[i] = area
            start_angle = end_angle
    return weights;


radial_bins = np.linspace(0,1,16)

#radius = input("enter the radius: ")
radius = np.float(sys.argv[1])
num_iter = 100000
proj_rad = np.zeros(num_iter,dtype='f4')
proj_x  = np.zeros_like(proj_rad)
proj_y  = np.zeros_like(proj_rad)
for i in range(0,num_iter):
    x = gauss(0,1)
    y = gauss(0,1)
    z = gauss(0,1)
    rr = x**2 + y**2 + z**2
    
    proj_rad[i]=radius*np.sqrt(x**2/rr + y**2/rr)
    #r = np.sqrt(rr)    
    #proj_x[i] = x/r
    #proj_y[i] = y/r


#plt.figure()
#plt.plot(proj_x,proj_y,'o',alpha=0.2)
#plt.axis('equal')
#plt.show()

H,x_bins = np.histogram(proj_rad,bins=radial_bins)
w = radial_weights(radius,radial_bins)
for i in range(0,radial_bins.size-1):
    #print i,radial_bins[i],"->",radial_bins[i+1]
    print '\t',i,"--",float(H[i])/float(num_iter), w[i]
print np.sum(H)/float(num_iter)
