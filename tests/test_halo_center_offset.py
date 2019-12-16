#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt


# def get_random_unit_vector():
#     xyz = np.random.normal(size=3)
#     xyz2 = xyz*xyz
#     r = np.sqrt(np.sum(xyz2))
#     return xyz/r

def get_random_unit_vectors(size):
    xyz = np.random.normal(size=(3,size))
    xyz2 = np.square(xyz)
    r = np.sqrt(np.sum(xyz2, axis=0))
    res = xyz/r
    return res

def get_radius_displacement(size, mean=45, width=15):
    print(size)
    rad = np.random.normal(loc=mean, scale=width, size=size)
    # if any of the radii are less than zero, replace them with new
    # ones that will be above zero via recursion. If a high fraction
    # of are below zero, this isn't a really an efficient solution.
    slct = rad<0
    if np.sum(slct)>0:
        rad[slct] = get_radius_displacement(np.sum(slct), mean=mean, width=width)
    return rad


def get_total_miscentering(size):
    xyz = get_random_unit_vectors(size)
    rad = get_radius_displacement(size)
    xyz_rad = xyz*rad
    return xyz_rad

def test_halo_center_offset():
    # xyzs = get_random_unit_vectors(10000000)
    # plt.figure()
    # plt.plot(xyzs[0,:], xyzs[1,:], ',b', alpha=0.002)
    # plt.axis('equal')
    # plt.figure()
    # plt.plot(xyzs[1,:], xyzs[2,:], ',b', alpha=0.002)
    # plt.axis('equal')
    # plt.figure()
    # plt.plot(xyzs[0,:], xyzs[2,:], ',b', alpha=0.002)
    # plt.axis('equal')

    # plt.show()
    # displacement = get_displacement(10000)
    # plt.figure()
    # h, xbins = np.histogram(displacement, bins=np.linspace(-100, 200, 100))
    # xbins_cen = (xbins[:-1]+xbins[1:])/2.0
    # plt.plot(xbins_cen, h)

    # plt.show()
    xyzs = get_total_miscentering(100000)
    plt.figure()
    plt.plot(xyzs[0,:], xyzs[1,:], ',b', alpha=0.2)
    plt.axis('equal')
    

if __name__ == "__main__":
    # get_random_unit_vectors(4)
    
    test_halo_center_offset()
    plt.show()
