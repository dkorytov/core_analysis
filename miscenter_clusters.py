#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

import dtk


def get_random_unit_vectors(size):
    xyz = np.random.normal(size=(3,size))
    xyz2 = np.square(xyz)
    r = np.sqrt(np.sum(xyz2, axis=0))
    res = xyz/r
    return res

def get_radius_displacement(size, mean, width):
    rad = np.random.normal(loc=mean, scale=width, size=size)
    # if any of the radii are less than zero, replace them with new
    # ones that will be above zero via recursion. If a high fraction
    # of are below zero, this isn't a really an efficient solution.
    slct = rad<0
    if np.sum(slct)>0:
        rad[slct] = get_radius_displacement(np.sum(slct), mean=mean, width=width)
    return rad

def get_total_miscentering(size, mean, width):
    xyz = get_random_unit_vectors(size)
    rad = get_radius_displacement(size, mean, width)
    xyz_rad = xyz*rad
    return xyz_rad

def load_clusters(fname):
    cat = {}
    hfile = h5py.File(fname, 'r')
    for key in hfile.keys():
        cat[key] = hfile[key][()]
    return cat

def save_clusters(fname, cat):
    print('saving...', fname)
    hfile = h5py.File(fname, 'w')
    for key in cat.keys():
        hfile[key] = cat[key]
        # print(key, cat[key])
    hfile.close()

def apply_miscentering_to_cluster(cat, miscenter_mean, miscenter_width):
    size = cat['fof_halo_center_x'].size
    xyz  = get_total_miscentering(size, miscenter_mean, miscenter_width)
    xx =  cat['fof_halo_center_x']
    # print(cat['fof_halo_center_x'])
    cat['fof_halo_center_x'] = cat['fof_halo_center_x']+xyz[0,:]
    # print(cat['fof_halo_center_x'])
    # print(cat['fof_halo_center_x']-xx)
    cat['fof_halo_center_y'] = cat['fof_halo_center_y']+xyz[1,:]
    cat['fof_halo_center_z'] = cat['fof_halo_center_z']+xyz[2,:]
     
def miscenter_clusters(param_fname):
    param = dtk.Param(param_fname)
    sod_input = param.get_string('sod_input')
    sod_output = param.get_string('sod_output')
    step = param.get_int('step')
    # re-scale from kpc to mpc
    miscenter_mean = param.get_float('miscenter_mean')/1000
    miscenter_width = param.get_float('miscenter_width')/1000
    sod_input = sod_input.replace('${step}', str(step))
    sod_output = sod_output.replace('${step}', str(step))
    assert sod_input != sod_output, 'Output file cannot be the same as the input file'
    cat = load_clusters(sod_input)
    apply_miscentering_to_cluster(cat, miscenter_mean, miscenter_width)
    save_clusters(sod_output, cat)
    

if __name__ == "__main__":
    miscenter_clusters(sys.argv[1])
