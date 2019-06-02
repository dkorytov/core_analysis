#!/usr/bin/env python2.7

from __future__ import print_function, division 

import numpy as np
import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys
import time
import h5py
from astropy.cosmology import WMAP7


def load_abundance_line(fname):
    print(fname)
    hfile = h5py.File(fname, 'r')
    print(hfile.keys())
    m_infall = hfile['abund_infall_mass'].value
    radius   = hfile['abund_radius'].value
    return m_infall, radius

def get_fname(param):
    core_loc = param.get_string('core_loc')
    step = param.get_int('step')
    expected_abundance = param.get_float("expected_comov_abundance")
    fname = "tmp_hdf5/{}/abundance={}.hdf5".format(core_loc, expected_abundance)
    return fname

def compare_abundance_line(param_fname1, param_fname2):
    param1 = dtk.Param(param_fname1)
    param2 = dtk.Param(param_fname2)
    m_infall1, radius1 = load_abundance_line(get_fname(param1))
    m_infall2, radius2 = load_abundance_line(get_fname(param2))
    f, ax = plt.subplots(1)
    ax.loglog(m_infall1, radius1, label='arg1')
    ax.loglog(m_infall2, radius2, label='arg2')
    plt.show()


if __name__ == "__main__":
    compare_abundance_line(sys.argv[1], sys.argv[2])
