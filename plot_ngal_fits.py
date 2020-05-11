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
    if file_name not in load_clusters._cache:
        cluster_data = ClusterData()
        cluster_data.load_file(file_name)
    else:
        cluster_data = load_clusters._cache[file_name]
    return cluster_data

load_clusters._cache = {}

def get_ngal_fit(param_fname, cluster_num, color, plot_fit=True, spider=False, manual_calc=False):
    param = dtk.Param(param_fname)
    cluster_loc = param.get_string('cluster_loc')
    if cluster_num is None:
        cluster_num = param.get_int('cluster_load_num')
    zmrh5_loc = param.get_string('zmrh5_loc')
    zmr_sdss  = ZMR(zmrh5_loc)
    zmr_fit = ZMR("output/"+param_fname+"/zmr_lkhd_cores.param")
    m_bins = zmr_fit.m_bins
    r_bins = zmr_fit.r_bins
    zmr_core_ngal, zmr_core_ngal_err = zmr_fit.get_ngal() # only one z-bin, so we don't select it out
    zmr_core_ngal = zmr_core_ngal[0]
    zmr_core_ngal_err = zmr_core_ngal_err[0]
    zmr_sdss_ngal, zmr_sdss_ngal_err = zmr_sdss.get_ngal()
    zmr_sdss_ngal = zmr_sdss_ngal[0]
    zmr_sdss_ngal_err = zmr_sdss_ngal_err[0]

    if manual_calc:
        model_fit_fname = "figs/"+param_fname+"/calc_likelihood_bounds.py/grid_fit_param.txt"
        model_fit = load_fit_limits(model_fit_fname)
        m_infall = 10**model_fit['mi']
        if 'rd' in model_fit:
            # print(model_fit['rd'])
            r_disrupt = model_fit['rd']/1000.0 #convert to mpc/h from kpc/h
        else:
            r_disrupt = np.inf
        # print("\ncalculating ngal for ", param_fname)
        # print("\tmodel_fit_fname:", model_fit_fname)
        # print("\tmodel params: {:.2e} {:.3f}".format(m_infall, r_disrupt))
        print(cluster_loc)
        cluster_data = load_clusters(cluster_loc)
        if cluster_num == -1:
            cluster_num = cluster_data.num
        cluster_ngal = np.zeros(cluster_num)
        cluster_m_i = np.zeros(cluster_num)
        for i in range(0, cluster_num):
            mass_index = cluster_data.get_cluster_mass_bin(i, m_bins)
            cluster_m_i[i] = mass_index
            cluster_ngal[i] = cluster_data.get_ngal(i, m_infall, r_disrupt)[1]
        ngal_mean = np.zeros(len(m_bins)-1)
        ngal_err = np.zeros(len(m_bins)-1)
        ngal_std = np.zeros(len(m_bins)-1)
        for i in range(0, len(m_bins)-1):
            slct = cluster_m_i == i
            ngal_mean[i] = np.mean(cluster_ngal[slct])
            ngal_std[i]  = np.std(cluster_ngal[slct])
            ngal_err[i]  = ngal_std[i]/np.sqrt(np.sum(slct))
            # print("{:.2e}->{:.2e}: {}".format(m_bins[i], m_bins[i+1], np.sum(slct)))
        plt.plot(dtk.bins_avg(m_bins), ngal_mean, '-x', color=color, label='Ngal recalc')
    if plot_fit:
        plt.plot(dtk.bins_avg(m_bins), zmr_core_ngal, '-', color=color)
        plt.fill_between(dtk.bins_avg(m_bins), zmr_core_ngal-zmr_core_ngal_err, zmr_core_ngal+zmr_core_ngal_err, color=color, alpha=0.3)
    offset_amount = 1.025
    if spider:
        markerfacecolor='None'
        markeredgecolor=color
        xaxis_offset=offset_amount
        lw = 1
    else:
        markerfacecolor=color
        markeredgecolor='None'
        xaxis_offset=1./offset_amount
        lw = 2
    
    # remove problematic 2.5 L* low mass cluster in the spider sample
    if "mstar-1" in param_fname and "spider" in param_fname:
        print("SPIDERSS!: ", zmr_sdss_ngal)
        zmr_sdss_ngal[zmr_sdss_ngal < 0.1 ] = np.nan
    plt.errorbar(dtk.bins_avg(m_bins)*xaxis_offset, zmr_sdss_ngal,
                 yerr=zmr_sdss_ngal_err, fmt='o', capsize=0, lw=lw, color=color,
                 markeredgecolor=markeredgecolor, markerfacecolor=markerfacecolor)
    # plt.fill_between(dtk.bins_avg(m_bins), ngal_mean-ngal_err, ngal_mean+ngal_err, color=color, alpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    # plt.legend(loc='best')
def format_plot():

    p4 = plt.plot([],[], 'tab:purple',  lw=5, label=r'{:1.2f}~L$_*$'.format(0.4))
    p3 = plt.plot([],[], 'tab:red', lw=5, label=r'{:1.2f}~L$_*$'.format(0.63))
    p2 = plt.plot([],[], 'tab:green', lw=5, label=r'{:1.2f}~L$_*$'.format(1.0))
    p12 = plt.plot([],[], 'tab:orange',lw=5, label=r'{:1.2f}~L$_*$'.format(1.58))
    p1 = plt.plot([],[],  'tab:blue',lw=5, label=r'{:1.2f}~L$_*$'.format(2.5))
    plt.errorbar([], [], yerr=[], fmt='o', lw=2, color='k', label="redMaPPer", capsize=0)
    plt.plot([], [], color='k', label="Core Model")
    # plt.errorbar([], [], yerr=[], fmt='o', lw=1, color='k', markerfacecolor='none', label='SPIDERS clusters', capsize=0)
    plt.legend(ncol=2, loc='best', framealpha=0.0)

    plt.xlabel(r'M$_{200c}$ [h$^{-1}$ M$_\odot$]')
    plt.ylabel(r'Projected N$_{\rm{gal}}$')
    plt.ylim([1e-1, 3e3])
    plt.xlim([1e14, 5e15])
    plt.tight_layout()

def plot_ngal_fits():
    get_ngal_fit("params/cfn/simet/mstar1/mean/a3_rd.param", None, 'c')
    get_ngal_fit("params/cfn/simet/mstar0.5/mean/a3_rd.param", None, 'g')
    get_ngal_fit("params/cfn/simet/mstar0/mean/a3_rd.param", None, 'b')
    get_ngal_fit("params/cfn/simet/mstar-1/mean/a3_rd.param", None, 'r')

    #just spider points
    get_ngal_fit("params/cfn/spider/mstar1/mean/spider_rd.param", None, 'c', plot_fit=False, spider=True)
    get_ngal_fit("params/cfn/spider/mstar0.5/mean/spider_rd.param", None, 'g', plot_fit=False, spider=True)
    get_ngal_fit("params/cfn/spider/mstar0/mean/spider_rd.param", None, 'b', plot_fit=False, spider=True)
    get_ngal_fit("params/cfn/spider/mstar-1/mean/spider_rd.param", None, 'r', plot_fit=False, spider=True)

    # get_ngal_fit("params/cfn/spider/mstar0/mean/spider_rd.param", None, 'm', plot_fit=False, spider=True)
    # get_ngal_fit("params/cfn/spider/mstar0/mean/bcg_rd.param", None, 'c', plot_fit=False, spider=True)
    format_plot()
def plot_ngal_fits2(pattern, mstars):
    color_cycle = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    for mstar, color in zip(mstars, color_cycle):
        get_ngal_fit(pattern.replace("${mstarval}", mstar), None, color)
    format_plot()

if __name__ == "__main__":
    if len(sys.argv) > 2:
        plot_name = sys.argv[1]
    else:
        plot_name = "OR_McClintock2019"
    mstars = ['-1', '-0.5', '0', '0.5', '1']
    if plot_name == "OR_Simet2017":
        pattern = 'params/rmba/auto/make_all_OR.high_richness.low_rez.min20.sh/crit/mstar${mstarval}/OR_rd_zoom.param'
        plot_ngal_fits2(pattern, mstars)
    elif plot_name == "OR_McClintock2019":
        pattern = 'params/rmba/auto/make_all_OR.McClintock.high_richness.low_rez.min20.sh/crit/mstar${mstarval}/OR_rd_zoom.param'
        plot_ngal_fits2(pattern, mstars)
    # plot_ngal_fits()
    dtk.save_figs("figs/"+__file__+"/"+plot_name+"/", extension='.pdf')
    plt.show()
