#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import os
import sys
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr

import dtk
import h5py
from matplotlib.patches import Circle

def get_cores(param_filename):
    cat = {}
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string("cluster_loc")
    hfile = h5py.File(cluster_loc, 'r')
    cat['radius'] = hfile['cores/core_r'].value
    cat['mass']   = hfile['cores/core_m'].value
    return cat

def plot_core_dist(param_filename):
    cores = get_cores(param_filename)
    mass_bins = np.logspace(10, 15, 5)
    mass_bins_cen = dtk.bins_avg(mass_bins)
    r_bins = np.logspace(-4,0, 64)
    r_bins_cen = dtk.bins_avg(r_bins)
    colors = plt.cm.copper(np.linspace(0,1,len(mass_bins_cen)))
    fig = plt.figure()
    for i in range(0, len(mass_bins_cen)):
        slct = (cores['mass'] > mass_bins[i]) & (cores['mass'] < mass_bins[i+1])
        h, _ = np.histogram(cores['radius'][slct], bins=r_bins, density=True)
        label = r"${:.2e}<M_{{200}}<{:.2e}$".format(mass_bins[i], mass_bins[i+1])
        plt.semilogx(r_bins_cen, h, color=colors[i], label=label)
    if "qc" in param_filename:
        plt.title("QContinuum Cores")
    else:
        plt.title("Outer Rim Cores")
    plt.ylabel("PDF")
    plt.xlabel("Core Radius [Mpc/h]")

    cmap = plt.cm.copper
    norm = clr.LogNorm(vmin=mass_bins_cen[0], vmax=mass_bins_cen[-1])
    cax, _ = matplotlib.colorbar.make_axes(plt.gca())
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb.set_label('Infall Mass [Msun/h]')
    plt.show()

if __name__ == "__main__":
    plot_core_dist(sys.argv[1])
