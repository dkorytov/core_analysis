#!/usr/bin/env python2.7

from __future__ import print_function, division
import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

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
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple']
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

    #offset for sdss vs spiders
    if spider_zmr is None: 
        offset = 0
    else:
        offset = .005
    plt.errorbar(r_bins_cen-offset, sdss_gal, yerr=+sdss_gal_err, fmt='o',
                 label='redMaPPer clusters', lw=2, color='k', capsize=0,
                 )
    if spider_zmr is not None:
        spider_gal = spider_zmr.zmr_gal_density[0, mass_i, :]
        spider_gal_err = spider_zmr.zmr_gal_density_err[0, mass_i, :]
        plt.errorbar(r_bins_cen+offset, spider_gal, yerr=+spider_gal_err, fmt='o',
                     label='SPIDERS clusters', lw=1, color='k', capsize=0, markerfacecolor='None')
       
    plt.yscale('log')
    plt.legend(loc='best', framealpha=0.0, ncol=2)
    title = r"{:.2f}$<$log$_{{10}}$(M$_{{200m}}$)$<${:.2f}".format(np.log10(sdss_zmr.m_bins[mass_i]), np.log10(sdss_zmr.m_bins[mass_i+1]))
    # plt.title(title)
    plt.text(0.05, 0.05, title, transform=plt.gca().transAxes)
    plt.xlabel(r'r/R$_{200}$')
    plt.ylabel(r"Galaxy Surface Density")
    plt.tight_layout()

def plot_profiles_mstar(param_base, mstar):
    # param_base = "params/cfn/simet/mstar0/mean/a3_@model@.param"
    # param_base = "params/rmba/auto/make_all_OR.McClintock.high_richness.low_rez.min20.sh/crit/mstar0/OR_@model@.param"
    param_base = param_base.replace("@mstar@", mstar)
    models = ["mi", "rd", "rm", "rd_rm",]
    model_names = ["Mi", "Rd", "Rm", "RdRm",]
    param_fnames = [param_base.replace("@model@", model) for model in models]
    model_zmrs = [load_zmr(param_fname) for param_fname in param_fnames]
    sdss_zmr = load_zmr(param_fnames[0], sdss=True) 
    # spider_zmr = load_zmr("params/cfn/spider/mstar0/mean/bcg.xray_mass.rd.param", sdss=True)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 0,)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 1,)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 2,)
    plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 3,)
    # plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 4,)
    # plot_multiple_model_profiles_one_mass(sdss_zmr, model_zmrs, model_names, 5,)

def plot_multiple_model_profiles(param_base, fig_dir):
    mstars = ['-1', '-0.5', '0', '0.5', '1']
    for mstar in mstars:
        plot_profiles_mstar(param_base, mstar)
        dtk.save_figs(fig_dir+"/mstar"+mstar+"/", extension=".pdf", reset_count=True)
        dtk.save_figs(fig_dir+"/mstar"+mstar+"/", extension=".png")
        plt.show()
        plt.close('all')

if __name__ == "__main__":
    param_bases = {"OR_McClintock":"params/rmba/auto/make_all_OR.McClintock.high_richness.low_rez.min20.sh/crit/mstar@mstar@/OR_@model@.param",
                    "OR_Simet":     "params/rmba/auto/make_all_OR.high_richness.low_rez.min20.sh/crit/mstar@mstar@/OR_@model@.param",
                    "OR_Baxter":    "params/rmba/auto/make_all_OR.Baxter.high_richness.low_rez.min20.sh/crit/mstar@mstar@/OR_@model@.param",
                    "OR_Farahi":    "params/rmba/auto/make_all_OR.Farahi.high_richness.low_rez.min20.sh/crit/mstar@mstar@/OR_@model@.param"}

    if len(sys.argv) >= 2:
        plot_names = syst.argv[1:]
    else:
        plot_names = list(param_bases.keys())
        
    for plot_name in plot_names:        
        plot_multiple_model_profiles(param_bases[plot_name], "figs/"+__file__+"/"+plot_name+"/")    
    plt.show()
