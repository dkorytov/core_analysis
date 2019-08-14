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
from core_fit2_util import *
from calc_ngal import *


n2merger= None

def test_core_catalog(core_loc, rank_num):
    print(core_loc)
    dtk.gio_inspect(core_loc)
    print("loading host htag")
    host_htag = dtk.gio_read(core_loc, 'fof_halo_tag', rank_num)
    print("loading infall step")
    infall_step  = dtk.gio_read(core_loc, 'infall_step', rank_num)
    print("loading infall htag")
    infall_htag  = dtk.gio_read(core_loc, 'infall_fof_halo_tag', rank_num)
    print("loading infall mass")
    infall_mass  = dtk.gio_read(core_loc, 'infall_mass', rank_num)
    print("loading central")
    core_central  = dtk.gio_read(core_loc, 'central', rank_num)
    
    mi = np.argmax(infall_mass)
    print(infall_mass[mi], np.log10(infall_mass[mi]))
    htag = host_htag[mi]
    slct = host_htag == htag
    host_htag = host_htag[slct]
    infall_htag = infall_htag[slct]
    infall_step = infall_step[slct]
    infall_mass = infall_mass[slct]
    core_central = core_central[slct]
    for i in range(0, np.sum(slct)):
        if infall_step[i] == 401:
            print(host_htag[i], infall_step[i], core_central[i])
    
def plot_saved_clusters(param_filename):
    global n2merger
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string("cluster_loc")
    cluster_data = ClusterData();
    n2lib_loc = "lib/libn2merg.so"
    n2merger = N2Merger(n2lib_loc)
    core_loc = param.get_string('core_loc').replace("${step}", str(401)).replace("_500L", "")
    fit_mi, fit_rd = get_fit_param(param_filename)
    print("Infall mass: {:.2e}\nR disrupt: {:.4f}".format(fit_mi, fit_rd))
    # test_core_catalog(core_loc, 1)
    

    print(n2lib_loc)
    cluster_data.load_file(cluster_loc)
    for i in range(0,cluster_data.num):
        if cluster_data.mass[i] > 1e15:
        # cluster_data.plot_fit_cluster(i, 12.2, 0.02)
            # cluster_data.plot_cluster(i)
            cluster_data.plot_fit_cluster(i, fit_mi, fit_rd)
            plt.show()

def plot_hist(data, bins, style='-', label=None, ):
    h, xbins = np.histogram(data, bins=bins)
    bins_cen = dtk.bins_avg(xbins)
    plt.plot(bins_cen, h, style, label=label)

def test_saved_clusters(param_filename):
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string('cluster_loc')
    sod_loc = param.get_string("sod_loc")
    sod_hdf5 = param.get_string("sod_hdf5")
    
    sod_data_raw = SODData();

    sod_data_saved = SODData();
    sod_data_saved.load_sod('/media/luna1/dkorytov/data/OuterRim/sod_200m/sod_m200m.401.test1.hdf5', True)

    saved_cluster = ClusterData()
    saved_cluster.load_file('tmp_hdf5/clusters_OR_M200m.test1.hdf5')
    # sod_data_raw.load_sod('/media/luna1/dkorytov/data/OuterRim/sod/m000.401.sodproperties', False)
    plt.figure()
    bins = np.logspace(12, 16, 32)
    plt.title("Mass")
    # plot_hist(sod_data_raw.mass, bins, label='sod_raw')
    plot_hist(sod_data_saved.mass, bins, label='sod saved')
    plot_hist(saved_cluster.mass, bins, label='saved cluster')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('m200')
    plt.legend(loc='best')

    plt.figure()
    bins = np.linspace(0,5,64)
    plt.title("radius")
    # plot_hist(sod_data_raw.rad, bins, label='sod_raw')
    plot_hist(sod_data_saved.rad, bins, label='sod saved')
    plot_hist(saved_cluster.rad, bins, label='saved cluster')
    plt.yscale('log')
    plt.xlabel('r200')
    plt.legend(loc='best')
    plt.show()

def plot_mass_function(param_filename):
    cluster_data = ClusterData(param_filename)
    mass_bins = np.logspace(14,16,16)
    mass_bins_cen = dtk.bins_avg(mass_bins)
    h, _ = np.histogram(clusters.mass, bins=mass_bins)
    plt.figure()
    plt.plot(mass_bins_cen, h)

    

if __name__ == "__main__":
    # plot_mass_function(sys.argv[1])
    plot_saved_clusters(sys.argv[1])
    #test_saved_clusters(sys.argv[1])


