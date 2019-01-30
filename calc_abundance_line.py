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
from core_fit2_util import *
from catalog_reader import Catalog,frag_to_real
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm


def load_core_cat(core_loc):
    print("loading: ", core_loc)
    t1 = time.time()
    dtk.gio_inspect(core_loc)
    cat = {}
    cat['infall_mass'] = dtk.gio_read(core_loc, 'infall_mass')
    cat['radius']      = dtk.gio_read(core_loc, 'radius')
    print(t1-time.time())
    return cat

def calc_abundance_infall_mass(cat, radius_cut, expected_number):
    infalls = cat['infall_mass'][cat['radius']<radius_cut]
    if len(infalls) < expected_number:
        return None
    else:
        srt = np.argsort(infalls)
        return infalls[srt][-expected_number]

def calc_abundance_radius(cat, infall_cut, expected_number):
    radius = cat['radius'][cat['infall_mass']>infall_cut]
    if len(radius) < expected_number:
        return None
    else:
        srt = np.argsort(radius)
        return radius[srt][expected_number]

def calc_abundance_line(param_fname, plot = False):
    param = dtk.Param(param_fname)
    core_loc = param.get_string('core_loc')
    rL = param.get_float('rL')
    expected_comov_abundance = param.get_float('expected_comov_abundance')
    step = param.get_int('step')
    cat = load_core_cat(core_loc.replace('${step}', str(step)))
    expected_num = int(rL*rL*rL*expected_comov_abundance)
    x = np.logspace(11,12.6,50)
    y = np.zeros_like(x)
    for i in range(0, len(x)):
        y[i] = calc_abundance_radius(cat, x[i], expected_num)
    if plot:
        print('num/vol: ', expected_comov_abundance )
        print('vol: ', rL*rL*rL)
        print('expected_num:', expected_num)
        
        print(x)
        print(y)
        plt.figure()
        plt.plot(x, y, '-x')
        plt.yscale('log')
        plt.xscale('log')
    y1 = np.logspace(-3, -.2, 25)
    x1 = np.zeros_like(y1)
    for i in range(0, len(x1)):
        x1[i] = calc_abundance_infall_mass(cat, y1[i], expected_num)
    if plot:
        plt.plot(x1, y1, '-x')
        print(x)
        print(y)

    abund_infall_mass = np.concatenate((x,x1))
    abund_radius = np.concatenate((y,y1))
    srt = np.argsort(abund_infall_mass)
    slct = np.isfinite(abund_infall_mass[srt]) & np.isfinite(abund_radius[srt])
    out_fname = 'tmp_hdf5/abundance={}.hdf5'.format(expected_comov_abundance)
    hfile = h5py.File(out_fname, 'w')
    hfile['abund_infall_mass'] = abund_infall_mass[srt][slct]
    hfile['abund_radius']=abund_radius[srt][slct]
    hfile.close()
    plt.show()
    


if __name__ == "__main__":
    calc_abundance_line(sys.argv[1], plot=True);
