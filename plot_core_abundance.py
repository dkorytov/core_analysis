#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr

import h5py
import sys

import dtk


def load_cores(core_fname, step):
    core_fname = core_fname.replace("${step}", str(step))
    dtk.gio_inspect(core_fname)
    result = {
        # "x": dtk.gio_read(core_fname, "x"),
        # "y": dtk.gio_read(core_fname, "y"),
        # "z": dtk.gio_read(core_fname, "z"),
        # "r": dtk.gio_read(core_fname, "r"),
        "infall_mass": dtk.gio_read(core_fname, "infall_mass"),
        "radius": dtk.gio_read(core_fname, "radius"),
    }
    return result

def plot_core_abundance(param_fname):
    param = dtk.Param(param_fname)
    cores = load_cores(param.get_string("core_loc"), param.get_int("step"))
    core_num = len(cores['infall_mass'])
    h, xbins, ybins = np.histogram2d(np.log10(cores['infall_mass']), np.log10(cores['radius']), bins =128)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap="Blues", norm=clr.LogNorm())
    plt.ylabel("Log10 Core Radius")
    plt.xlabel("Log10 Infall Mass")



    h_sum = np.cumsum(h, axis=1)
    # h_sum = h_sum[-1,:] - h_sum
    # h_sum = h_sum.T
    plt.figure()
    plt.pcolor(xbins, ybins, h_sum.T, cmap="Blues", norm=clr.LogNorm())
    plt.ylabel("Log10 Core Radius")
    plt.xlabel("Log10 Infall Mass")


    h_sum2 = np.cumsum(h_sum, axis=0)
    h_sum2 = h_sum2[-1,:] - h_sum2
    plt.figure()
    plt.pcolor(xbins, ybins, h_sum2.T, cmap="Blues", norm=clr.LogNorm())
    plt.contour(dtk.bins_avg(xbins), dtk.bins_avg(ybins), h_sum2.T)
    plt.ylabel("Log10 Core Radius")
    plt.xlabel("Log10 Infall Mass")

    plt.figure()
    plt.pcolor(xbins, ybins, h_sum2.T, cmap="Blues", norm=clr.LogNorm())
    plt.contour(dtk.bins_avg(xbins), dtk.bins_avg(ybins), h_sum2.T/(500*500*500), levels = [0.0010, 0.0015, 0.002])
    plt.ylabel("Log10 Core Radius")
    plt.xlabel("Log10 Infall Mass")

    plt.figure()
    percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    result1 = dtk.binned_percentile(np.log10(cores['infall_mass']), np.log10(cores['radius']), xbins, percentiles, minimum_number=100)
    print(result1.shape)
    plt.pcolor(xbins, ybins, h.T, cmap="Blues", norm=clr.LogNorm())
    plt.colorbar(label='Population Density')
    plt.plot(12.25, np.log10(0.05), "*", label ='model fit')
    for i, percentile in enumerate(percentiles):
        plt.plot(dtk.bins_avg(xbins), result1[:,i], '-', label="{}%".format(percentile))

    plt.legend(loc='best', framealpha=0)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        plot_core_abundance("params/cfn/simet/simet_rd_rm_wide_rdlock.param")
    else:
        plot_core_abundance(sys.argv[1])
