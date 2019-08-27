#!/usr/bin/env python2.7

from __future__ import print_function, division

from matplotlib import rc
# rc('font',**{'family': 'serif',
#              'serif':  ['DejaVu'],
#              'size':   15})
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np

from zmr import ZMR


def load_zmr(param_fname, sdss=False):
    if not sdss:
        fname = "output/"+param_fname+"/zmr_lkhd_cores.param"
    if sdss:
        fname = "output/"+param_fname+"/zmr_sdss.param"
    return ZMR(fname)

    
def plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, mass_i, spider_zmr=None):
    plt.figure()
    model_size = len(model_zmrs)
    colors = ['r', 'b', 'g', 'c']
    for i in range(0, model_size):
        gal_density = model_zmrs[i].zmr_gal_density[0, mass_i, :]
        gal_density_err = model_zmrs[i].zmr_gal_density_err[0, mass_i, :]
        r_bins_cen = dtk.bins_avg(model_zmrs[i].r_bins)
        plt.plot(r_bins_cen, gal_density, label=model_names[i], color=colors[i], lw=2)
        plt.fill_between(r_bins_cen, gal_density-gal_density_err, gal_density+gal_density_err, color=colors[i], alpha=0.2)

    sdss_gal = sdss_zmr.zmr_gal_density[0, mass_i, :]
    sdss_gal_err = sdss_zmr.zmr_gal_density_err[0, mass_i, :]
    # plt.plot(r_bins_cen, sdss_zmr.zmr_gal_density[0, mass_i, :], 'k', label='SDSS', lw=2)
    # plt.fill_between(r_bins_cen, sdss_gal-sdss_gal_err, sdss_gal+sdss_gal_err, color='k', alpha=0.2)
    plt.errorbar(r_bins_cen, sdss_gal, yerr=+sdss_gal_err, fmt='o',
                 label='redMaPPer Galaxies', lw=2, color='k', capsize=0,
                 )
    if spider_zmr is not None:
        spider_gal = spider_zmr.zmr_gal_density[0, mass_i, :]
        spider_gal_err = spider_zmr.zmr_gal_density_err[0, mass_i, :]
        plt.errorbar(r_bins_cen, spider_gal, yerr=+spider_gal_err, fmt='o',
                     label='SPIDER Galaxies', lw=2, color='k', capsize=0, markerfacecolor='None')
       
    plt.yscale('log')
    plt.legend(loc='best', framealpha=0.0, ncol=2)
    title = r"{:.2f}$<$log$_{{10}}$(M$_{{200m}}$)$<${:.2f}".format(np.log10(sdss_zmr.m_bins[mass_i]), np.log10(sdss_zmr.m_bins[mass_i+1]))
    # plt.title(title)
    plt.text(0.05, 0.05, title, transform=plt.gca().transAxes)
    plt.xlabel(r'r/R$_{200}$')
    plt.ylabel(r"Galaxy Surface Density")

def plot_mstar0():
    param_base = "params/cfn/simet/mstar0/mean/a3_@model@.param"
    models = ["mi", "rd", "rm", "rd_rm",]
    model_names = ["Mi", "Rd", "Rm", "RdRm",]
    param_fnames = [param_base.replace("@model@", model) for model in models]
    model_zmrs = [load_zmr(param_fname) for param_fname in param_fnames]
    sdss_zmr = load_zmr(param_fnames[0], sdss=True) 
    spider_zmr = load_zmr("params/cfn/simet/mstar0/mean/spider_rd.param", sdss=True)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 0, spider_zmr=spider_zmr)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 1, spider_zmr=spider_zmr)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 2, spider_zmr=spider_zmr)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 3, spider_zmr=spider_zmr)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 4, spider_zmr=spider_zmr)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 5, spider_zmr=spider_zmr)


def plot_multiple_model_profiles():
    plot_mstar0()


if __name__ == "__main__":
    plot_multiple_model_profiles()    
    dtk.save_figs("figs/"+__file__+"/", extension=".pdf")
    plt.show()
