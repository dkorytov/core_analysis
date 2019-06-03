#!/usr/bin/env python2.7

from __future__ import print_function, division 

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as  clr
import h5py 
import pandas as pd

import dtk

def load_cores(fname):
    print(fname)
    cat = {}
    hfile = h5py.File(fname, 'r')
    keys = ['m_infall', 'radius', 'm_peak', 'r_peak', 'r_peak_raw']
    for k in keys:
        cat[k] = hfile[k].value
    return cat


def check_fit_cores(fname, m_infall, r_radius):
    cat = load_cores(fname)
    slct1 = (cat['m_infall'] > m_infall) & (cat['radius'] < r_radius)
    slct2 = (cat['m_peak'] > m_infall) & (cat['radius'] < r_radius)
    slct3 = (cat['m_infall'] > m_infall) & (cat['r_peak'] < r_radius)
    slct4 = (cat['m_infall'] > m_infall) & (cat['r_peak_raw'] < r_radius)
    slct5 = (cat['m_peak'] > m_infall) & (cat['r_peak'] < r_radius)
    print("Default:   ", np.sum(slct1)/np.sum(slct1))
    print("m peak:    ", np.sum(slct2)/np.sum(slct1))
    print("r peak:    ", np.sum(slct3)/np.sum(slct1))
    print("r peak raw:", np.sum(slct4)/np.sum(slct1))
    print("r m peak:  ", np.sum(slct5)/np.sum(slct1))


if __name__ == "__main__":
    check_fit_cores('/media/luna1/dkorytov/data/OuterRim/cores_peak/OuterRim.401.central.hdf5', 10**12.2, 0.02)
