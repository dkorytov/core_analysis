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
from matplotlib.colors import LogNorm
import dtk
import sys
import time

from scipy.optimize import minimize
from plot_saved_clusters import ClusterData
from generate_parameter_dist import load_fit_limits

def get_fit_limits_fname(param_fname):
    result = "figs/{}/calc_likelihood_bounds.py/grid_fit_param.txt".format(param_fname)
    return result



def calc_ngal(param_fname):
    param = dtk.Param(param_fname)
    cluster_loc = param.get_string('cluster_loc')
    central = param.get_bool('force_central_as_galaxy')
    clusters = ClusterData();
    clusters.load_file(cluster_loc, treat_centrals = central, step=401)
    fit_fname = get_fit_limits_fname(param_fname)
    print(fit_fname)
    fit_param = load_fit_limits(fit_fname)
    fit_mi = 10**fit_param['mi']
    fit_rd = fit_param['rd']*1e-3
    print(fit_mi, fit_rd)
    ngals = []
    ngals_proj = []
    masses = []
    print(clusters.num)
    for i in range(0,clusters.num):
        ngal, ngal_proj = clusters.get_ngal(i, fit_mi, fit_rd)
        masses += [clusters.mass[i]]
        ngals += [ngal]
        ngals_proj+= [ngal_proj]
    masses = np.array(masses)
    ngals = np.array(ngals)
    ngals_proj = np.array(ngals_proj)
    xbins = np.logspace(14,15.5, 32)
    xbins_cen = dtk.bins_avg(xbins)
    plt.figure()
    plt.plot(masses, ngals, 'ob', alpha=0.3, label = 'Ngal')
    plt.plot(masses, ngals_proj, 'og', alpha=0.3, label='Projected Ngal')
    ngal_avg = dtk.binned_average(masses, ngals, xbins)
    ngal_proj_avg = dtk.binned_average(masses, ngals_proj, xbins)

    plt.plot(xbins_cen, ngal_avg, 'b', lw=2)
    plt.plot(xbins_cen, ngal_proj_avg, 'g', lw=2)

    plt.xscale('log')
    plt.show()

if __name__ == "__main__":
    calc_ngal(sys.argv[1])
