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
from scipy.special import erf


from scipy.optimize import minimize
from plot_saved_clusters import ClusterData
from generate_parameter_dist import load_fit_limits

class HODModel:
    def __init__(self, log_M_min, sigma, log_M_0, log_M_1, alpha):
        self.M_min = 10**log_M_min
        self.log_M_min = log_M_min
        self.sigma = sigma
        self.M_0 = 10**log_M_0
        self.M_1 = 10**log_M_1
        self.alpha = alpha
        # self.M_1a = 10**log_M_1a
        # self.b_ga = b_ga

    def hod_cen(self, M200m):
        return 0.5*(1.0+erf((np.log10(M200m) - self.log_M_min)/self.sigma))
        
    def hod_sat(self, M200m):
        result = self.hod_cen(M200m)*(((M200m - self.M_0)/self.M_1)**self.alpha)
        result[~np.isfinite(result)] = 0
        return result

    def hod(self, M200m):
        return self.hod_cen(M200m) + self.hod_sat(M200m)
                      

def get_zeng07_hod(Mr=None):
    if Mr is None:
        raise KeyError
    if Mr == -18:
        return HODModel(11.35, 0.25, 11.20, 12.40, 0.83)
    elif Mr == -18.5:
        return HODModel(11.46, 0.24, 10.59, 12.68, 0.97)
    elif Mr == -19:
        return HODModel(11.60, 0.26, 11.49, 12.83, 1.02)
    elif Mr == -19.5:
        return HODModel(11.75, 0.28, 11.69, 13.01, 1.06)
    elif Mr == -20.0:
        return HODModel(12.02, 0.26, 11.38, 13.31, 1.06)
    elif Mr == -20.5:
        return HODModel(12.30, 0.21, 11.84, 13.58, 1.12)
    elif Mr == -21.0:
        return HODModel(12.79, 0.39, 11.92, 13.94, 1.15)
    elif Mr == -21.5:
        return HODModel(13.38, 0.51, 13.94, 13.91, 1.04)
    elif Mr == -22.0:
        return HODModel(14.22, 0.77, 14.00, 14.69, 0.87)
    else:
        raise KeyError

def get_zeng07_Mrs():
    return [-18.0, -18.5, -19.0, -19.5, -20.0, -20.5, -21.0, -21.5, -22.0]

def get_fit_limits_fname(param_fname):
    result = "figs/{}/calc_likelihood_bounds.py/grid_fit_param.txt".format(param_fname)
    return result

def get_mr_limits(param_fname):
    if "/mstar0/" in param_fname:
        return -21.5, -21.0
    elif "/mstar0.5/" in param_fname:
        return -21.0, -20.5
    elif "/mstar1/" in param_fname:
        return -20.5, -20.0
    elif "/mstar-1/" in param_fname:
        return -22.5, -22.

def get_clusters(param_fname):
    param = dtk.Param(param_fname)
    cluster_loc = param.get_string('cluster_loc')
    central = param.get_bool('force_central_as_galaxy')
    if central:
        print("It's compact central")
    clusters = ClusterData();
    clusters.load_file(cluster_loc, treat_centrals = central, step=401)
    return clusters

def calc_ngal(param_fname):
    clusters = get_clusters(param_fname)
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
        ngal, ngal_proj = clusters.get_ngal(i, fit_mi, fit_rd, compact_central = central)
        masses += [clusters.mass[i]]
        ngals += [ngal]
        ngals_proj+= [ngal_proj]
    masses = np.array(masses)
    ngals = np.array(ngals)
    ngals_proj = np.array(ngals_proj)
    xbins = np.logspace(14,15.5, 32)
    xbins_cen = dtk.bins_avg(xbins)
    plt.figure()
    # plt.plot(masses, ngals, ',b', alpha=0.01, label = 'Ngal')
    # plt.plot(masses, ngals_proj, ',g', alpha=0.01, label='Projected Ngal')
    ngal_avg = dtk.binned_average(masses, ngals, xbins)
    ngal_proj_avg = dtk.binned_average(masses, ngals_proj, xbins)

    plt.plot(xbins_cen, ngal_avg, 'b', lw=2, label = "Core Ngal")
    plt.plot(xbins_cen, ngal_proj_avg, 'g', lw=2, label = "Core Ngal Avg Projected")
    plt.xscale('log')
    plt.yscale('log')

    ylim = plt.ylim()
    xlim = plt.xlim()
    mr_upper, mr_lower = get_mr_limits(param_fname)
    mass = np.logspace(11,15.5,100)
    model_lower = get_zeng07_hod(mr_lower)
    hod_tot_lower = model_lower.hod(mass)
    plt.plot(mass, hod_tot_lower, '--r', lw = 1, label="Mr<{:.1f} HOD".format(mr_lower))
    model_upper = get_zeng07_hod(mr_upper)
    hod_tot_upper = model_upper.hod(mass)
    plt.plot(mass, hod_tot_upper, 'r--', lw = 1, label="Mr<{:.1f} HOD".format(mr_upper))    
    plt.plot(mass, calc_log_mean(hod_tot_lower, hod_tot_upper), 'r', label='Mr<{:.1f} HOD'.format((mr_upper+mr_lower)/2.0))
    plt.legend(loc='best')

    plt.xlim(xlim)
    plt.ylim(ylim)

    
    plt.show()

def calc_zeng07():
    mass = np.logspace(11,15.5,100)
    Mrs = get_zeng07_Mrs()
    plt.figure()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#ff7f0e', '#8c564b']
    for i, Mr in enumerate(Mrs):
        hod_model = get_zeng07_hod(Mr)
        hod_cen = hod_model.hod_cen(mass)
        hod_sat = hod_model.hod_sat(mass)
        hod_tot = hod_model.hod(mass)
        plt.plot(mass, hod_tot, lw=2, color=colors[i])
        plt.plot(mass, hod_cen, lw=1, color=colors[i], alpha=0.3, ls='--')
        plt.plot(mass, hod_sat, lw=1, color=colors[i], alpha=0.3, ls=':')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

def calc_log_mean(data1, data2):
    a = np.sqrt(data1*data2)
    return a
def calc_

if __name__ == "__main__":
    calc_ngal(sys.argv[1])
    # calc_zeng07()
