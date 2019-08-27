#!/usr/bin/env python2.7

from __future__ import print_function, division

from matplotlib import rc
rc('font',**{'family': 'serif',
             'serif':  ['DejaVu'],
             'size':   15})
rc('text', usetex=True)

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np
import h5py

def test_infall_peak(fname):
    hfile = h5py.File(fname, 'r')
    m_infall = hfile['m_infall'].value
    m_peak   = hfile['m_peak'].value
    print(m_infall)
    print(m_peak)
    print(np.sum(m_infall == 0))
    print(np.sum(m_peak == 0))
    print(np.average(m_peak/m_infall))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        fname = '/media/luna1/dkorytov/data/OuterRim/cores_peak/OuterRim.401.central.hdf5'
    else:
        param = dtk.Param(sys.argv[1])
        fname = param.get_string(core_loc)
    test_infall_peak(fname)
