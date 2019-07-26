#!/usr/bin/env python2.7

from __future__ import print_function, division

from matplotlib import rc
rc('font',**{'family': 'serif',
             'serif':  ['Palatino'],
             'size':   15})
rc('text', usetex=True)

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

from zmr import ZMR
from generate_parameter_dist import load_fit_limits

# def load_zmr(param_fname, sdss=False):
#     if not sdss:
#         fname = "output/"+param_fname+"/zmr_lkhd_cores.param"
#     if sdss:
#         fname = "output/"+param_fname+"/zmr_sdss.param"
#     return ZMR(fname)


# def load_zmrs(param_fnames, sdss=False):
#     return [ load_zmr(param_fname, sdss) for param_fname in param_fnames]

def load_fit_limits_set(param_fnames):
    return [ load_fit_limits("figs/"+param_fname+"/calc_likelihood_bounds.py/grid_fit_param.txt") for param_fname in param_fnames]

    
def plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, mass_i):
    plt.figure()
    model_size = len(model_zmrs)
    colors = ['r', 'b', 'g', 'c']
    for i in range(0, model_size):
        gal_density = model_zmrs[i].zmr_gal_density[0, mass_i, :]
        gal_density_err = model_zmrs[i].zmr_gal_density_err[0, mass_i, :]
        r_bins_cen = dtk.bins_avg(model_zmrs[i].r_bins)
        plt.plot(r_bins_cen, gal_density, label=model_names[i], color=colors[i])
        plt.fill_between(r_bins_cen, gal_density-gal_density_err, gal_density+gal_density_err, color=colors[i], alpha=0.3)
    sdss_gal = sdss_zmr.zmr_gal_density[0, mass_i, :]
    sdss_gal_err = sdss_zmr.zmr_gal_density_err[0, mass_i, :]
    print(np.shape(sdss_gal))
    print(np.shape(sdss_gal_err))
    plt.plot(r_bins_cen, sdss_zmr.zmr_gal_density[0, mass_i, :], 'k', label='SDSS', lw=2)
    plt.fill_between(r_bins_cen, sdss_gal-sdss_gal_err, sdss_gal+sdss_gal_err, color='k', alpha=0.3)
    plt.yscale('log')
    plt.legend(loc='best')
    title = r"{:.2f}$<$log$_{{10}}$(M$_{{200m}}$)$<${:.2f}".format(np.log10(sdss_zmr.m_bins[i]), np.log10(sdss_zmr.m_bins[i+1]))
    # plt.title(title)
    plt.text(0.5, 0.9, title, horizontalalignment='center',
         verticalalignment='center', transform=plt.gca().transAxes)
    plt.xlabel(r'r/R$_{200}$')
    plt.ylabel(r"Galaxy Surface Density")



def plot_luminosity_dependence(param_base, mstars, model_params):
    # param_base = 

    # if "rm" in param_base:

    # else:
    #     model_params = ['mi', 'rd', 'x2']
    plt.figure()
    parameter_labels = {'mi': r"M$_{infall}$",
                        'rd': r"R$_{disrupt}$",
                        'rm': r"R$_{merge}$",
                        'x2': r"$\chi^2_{reduced}$",}
    param_fnames = [param_base.replace("@val@", str(mstar)) for mstar in mstars]
    # zmrs = load_zmrs(param_fnames)
    fits = load_fit_limits_set(param_fnames)
    gs = gridspec.GridSpec(len(model_params),1, hspace=0.1)
    share_x = None
    for i, model_param in enumerate(model_params):
        print(i)
        ax = plt.subplot(gs[i,0], sharex=share_x)
        model_param_vals = [fit[model_param] for fit in fits]
        ax.plot(mstars, model_param_vals)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.set_ylabel(parameter_labels[model_param])
        if model_param != "mi":
            ylim = ax.get_ylim()
            ax.set_ylim([0,ylim[1]])
        if i == len(model_params)-1:
            ax.set_xlabel("Galaxy Mstar Offset")
        else:
            ax.xaxis.set_ticklabels([])

def plot_multiple_model_profiles():
    plot_mstar0()

def init_luminosity_dependence_plot():
    plots_axies = {}
    gs = gridspec.GridSpec(4, 1, hspace=0.1)
    fig = plt.figure(figsize=(10,8))
    parameter_labels = {'mi': r"M$_{infall}$",
                        'rd': r"R$_{disrupt}$",
                        'rm': r"R$_{merge}$",
                        'x2': r"$\chi^2_{reduced}$",}
    shared_x = None
    ax_dict = {}
    for i, model_param in enumerate(['mi', 'rd', 'rm', 'x2']):
        ax = plt.subplot(gs[i,0], sharex=shared_x)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.set_ylabel(parameter_labels[model_param])
        if model_param != "mi":
            ylim = ax.get_ylim()
            ax.set_ylim([0,ylim[1]])
        if i == 3:
            ax.set_xlabel("Galaxy Mstar Offset")
        else:
            ax.xaxis.set_ticklabels([])
        ax_dict[model_param] = ax
    return fig, ax_dict
    

def plot_luminosity_dependence_single(fig, ax_dict, mstars, param_base, color, param_list):
    param_fnames = [param_base.replace("@val@", str(mstar)) for mstar in mstars]
    fits = load_fit_limits_set(param_fnames)
    for model_param in param_list:
        values = [fits[i][model_param] for i in range(0, len(mstars))]
        ax_dict[model_param].plot(mstars, values, color=color)
        
    
def plot_luminosity_dependence_all():
    fig, ax_dict = init_luminosity_dependence_plot()
    mstars = [-1, 0, 0.5, 1]
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a2_mi.param", 'r', ['mi', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a2_rd.param", 'b', ['mi', 'rd', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a2_rm.param", 'g', ['mi', 'rm', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a2_rd_rm.param", 'c', ['mi', 'rd', 'rm', 'x2'])
    plt.show()

if __name__ == "__main__":
    mstars = [-1, 0, 0.5, 1]
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a2_mi.param", mstars, ['mi', 'rd', 'rm', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rd.param", mstars, ['mi', 'rd', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rm.param", mstars, ['mi', 'rm', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rd_rm.param", mstars, ['mi', 'rd', 'rm', 'x2'])
    plot_luminosity_dependence_all()
    dtk.save_figs("figs/"+__file__+"/", extension=".pdf")
    plt.show()
