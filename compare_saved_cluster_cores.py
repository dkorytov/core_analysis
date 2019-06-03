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
    keys = ['core_r', 'core_m']
    for k in keys:
        cat[k] = hfile['cores/'+k].value
    return cat

def get_passed(cat, m_cut, r_cut):
    slct = (cat['core_m']>m_cut) & (cat['core_r'] < r_cut)
    return np.sum(slct)
    
def compare_saved_cluster_cores(fname1, fname2, m_cut, r_cut):
    cat1 = load_cores(fname1)
    cat2 = load_cores(fname2)
    num1 = get_passed(cat1, m_cut, r_cut)
    num2 = get_passed(cat2, m_cut, r_cut)
    print(fname1, num1, num1/num1)
    print(fname2, num2, num2/num1)

if __name__ == "__main__":
    compare_saved_cluster_cores(sys.argv[1], sys.argv[2], 10**12.2, 0.02)
