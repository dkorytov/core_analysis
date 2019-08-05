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
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize

from calc_ngal import *
from generate_parameter_dist import *
from zmr import ZMR
from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

def load_clusters(file_name):
    print("we are loading clusters")
    file_name = "tmp_hdf5/clusters_OR.hdf5"
    if file_name not in load_clusters._cache:
        cluster_data = ClusterData()
        cluster_data.load_file(file_name)
    else:
        cluster_data = load_clusters._cache[file_name]
    return cluster_data

load_clusters._cache = {}

def get_ngal_fit(param_fname, cluster_num, color):
    param = dtk.Param(param_fname)
    zmr_fit = ZMR("output/"+param_fname+"/zmr_lkhd_cores.param")
    m_bins = zmr_fit.m_bins
    print(m_bins)
    r_bins = zmr_fit.r_bins
    model_fit_fname = "figs/"+param_fname+"/calc_likelihood_bounds.py/grid_fit_param.txt"
    model_fit = load_fit_limits(model_fit_fname)
    m_infall = 10**model_fit['mi']
    if 'rd' in model_fit:
        print(model_fit['rd'])
        r_disrupt = model_fit['rd']/1000.0 #convert to mpc/h from kpc/h
    else:
        r_disrupt = np.inf
    print("\ncalculating ngal for ", param_fname)
    print("\tmodel_fit_fname:", model_fit_fname)
    print("\tmodel params: {:.2e} {:.3f}".format(m_infall, r_disrupt))
    
    cluster_data = load_clusters(param.get_string('cluster_loc'))
    cluster_ngal = np.zeros(cluster_num)
    cluster_m_i = np.zeros(cluster_num)
    for i in range(0, cluster_num):
        mass_index = cluster_data.get_cluster_mass_bin(i, m_bins)
        cluster_m_i[i] = mass_index
        cluster_ngal[i] = cluster_data.get_ngal(i, m_infall, r_disrupt)[0]
    ngal_mean = np.zeros(len(m_bins)-1)
    ngal_err = np.zeros(len(m_bins)-1)
    ngal_std = np.zeros(len(m_bins)-1)
    for i in range(0, len(m_bins)-1):
        slct = cluster_m_i == i
        ngal_mean[i] = np.mean(cluster_ngal[slct])
        ngal_std[i]  = np.std(cluster_ngal[slct])
        ngal_err[i]  = ngal_std[i]/np.sqrt(np.sum(slct))
        print("{:.2e}->{:.2e}: {}".format(m_bins[i], m_bins[i+1], np.sum(slct)))
    plt.plot(dtk.bins_avg(m_bins), ngal_mean, '-x', color=color)
    plt.yscale('log')
    plt.xscale('log')


def plot_ngal_fits():
    get_ngal_fit("params/cfn/simet/mstar-1/mean/a3_rd.param", 5000, 'r')
    get_ngal_fit("params/cfn/simet/mstar0/mean/a3_rd.param", 5000, 'b')
    get_ngal_fit("params/cfn/simet/mstar0.5/mean/a3_rd.param", 5000, 'g')
    get_ngal_fit("params/cfn/simet/mstar1/mean/a3_rd.param", 5000, 'c')

if __name__ == "__main__":
    plot_ngal_fits()
    plt.show()
