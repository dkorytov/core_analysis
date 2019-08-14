#!/usr/bin/env python2.7

from __future__ import print_function, division


import dtk
import sys
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

from zmr import ZMR
from generate_parameter_dist import load_fit_limits
from plot_ngal_fits import load_clusters
from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 
              'serif':['Computer Modern Roman'], 
              'size':18})

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
    fig = plt.figure(figsize=(8,8))
    parameter_labels = {'mi': r"M$_{in fall}$"+"\n"+r"[h$^{-1}$M$_{\odot}$]",
                        'rd': r"R$_{disrupt}$"+"\n"+r"[h$^{-1}$kpc]",
                        'rm': r"R$_{merge}$"+"\n"+r"[h$^{-1}$Mpc]",
                        'x2': r"$\chi^2_{reduced}$",}
    shared_x = None
    ax_dict = {}
    for i, model_param in enumerate(['mi', 'rd', 'rm', 'x2']):
        ax = plt.subplot(gs[i,0], sharex=shared_x)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.set_ylabel(parameter_labels[model_param])
        if i == 3:
            ax.set_xlabel("Galaxy Mstar Offset")
        else:
            ax.xaxis.set_ticklabels([])
        if model_param == 'x2' or model_param == 'mi':
            ax.set_yscale('log')

        if i == 0:
            ax.plot([], [], 'r', lw=3, label='Mi')
        elif i == 1:
            ax.plot([], [], 'b', lw=3, label='Rd')
        elif i == 2:
            ax.plot([], [], 'g', lw=3, label='Rm')
        elif i == 3:
            ax.plot([], [], 'c', lw=3, label='RdRm')
        ax.legend(loc='best', framealpha=0)
        ax_dict[model_param] = ax
        
    return fig, ax_dict, gs
    
def plot_luminosity_dependence_single(fig, ax_dict, mstars, param_base, color, param_list):
    param_fnames = [param_base.replace("@val@", str(mstar)) for mstar in mstars]
    fits = load_fit_limits_set(param_fnames)
    [ print(param_fname, fit) for param_fname, fit in zip(param_fnames, fits) ] 
    for model_param in param_list:
        values = np.array([fits[i][model_param] for i in range(0, len(mstars))])
        if model_param == 'mi':
            ax_dict[model_param].plot(mstars, 10**(values), color=color)
        else:
            ax_dict[model_param].plot(mstars, values, color=color)
        if model_param is not "x2":
            values_upper_err = np.array([fits[i][model_param+"_upper_err"] for i in range(0, len(mstars))])
            values_lower_err = np.array([fits[i][model_param+"_lower_err"] for i in range(0, len(mstars))])
            if model_param == 'mi':
                ax_dict[model_param].fill_between(mstars, 10**(values+values_upper_err), 10**(values-values_lower_err), color=color, alpha=0.3)
            else:
                ax_dict[model_param].fill_between(mstars, values+values_upper_err, values-values_lower_err, color=color, alpha=0.3)

        print(model_param, values)

def plot_luminosity_dependence_parameters_all():
    fig, ax_dict, gs = init_luminosity_dependence_plot()
    mstars = [-1, 0, 0.5, 1]
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a3_mi.param", 'r', ['mi', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a3_rd.param", 'b', ['mi', 'rd', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/a3_rm.param", 'g', ['mi', 'rm', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars,
                                      "params/cfn/simet/mstar@val@/mean/a3_rd_rm.param", 'c', ['mi',
                                                                                               'rd', 'rm', 'x2'])
    gs.tight_layout(fig)
    
    fig, ax_dict, gs = init_luminosity_dependence_plot()
    mstars = [0, 0.5, 1]
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/crit/a3_mi.param", 'r', ['mi', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/crit/a3_rd.param", 'b', ['mi', 'rd', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/crit/a3_rm.param", 'g', ['mi', 'rm', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/crit/a3_rd_rm.param", 'c', ['mi', 'rd', 'rm', 'x2'])
    gs.tight_layout(fig)

    fig, ax_dict, gs = init_luminosity_dependence_plot()
    mstars = [0, 0.5, 1]
    # plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/qc_mi.param", 'r', ['mi', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/qc_rd.param", 'b', ['mi', 'rd', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/qc_rm.param", 'g', ['mi', 'rm', 'x2'])
    plot_luminosity_dependence_single(fig, ax_dict, mstars, "params/cfn/simet/mstar@val@/mean/qc_rd_rm.param", 'c', ['mi', 'rd', 'rm', 'x2'])
    gs.tight_layout(fig)


def plot_luminosity_dependent_ngal_all(param_file):
    param = dtk.Param(param_file)
    sdss_zmr_loc = param.get_string('zmrh5_loc')
    print(sdss_zmr_loc)
    hfile = h5py.File(sdss_zmr_loc, 'r')
    zmr_sdss = ZMR(file_loc="output/"+param_file+"/zmr_sdss.param")
    if dtk.file_exists("output/"+param_file+"/zmr_lkhd_cores.param"):
        fname = "output/"+param_file+"/zmr_lkhd_cores.param"
        zmr_cores = ZMR(file_loc=fname)
        print("likelihood zmrs")
    else:
        fname = "output/"+param_file+"/zmr_cores.param"
        zmr_cores = ZMR(file_loc=fname)
        print("fit zmrs")

    mass_bins = zmr_cores.m_bins 
    mass_bins_centers = dtk.bins_avg(mass_bins)
    zmr_Ngal = zmr_cores.zmr_gal_counts[0].sum(axis=1)/zmr_cores.zm_counts[0]
    

def get_survival_rate_cluster(cluster_data, cluster_num, mass_bins, model_fit, ):
    m_infall = 10**model_fit['mi']
    r_disrupt = model_fit['rd']/1000
    if cluster_num == -1:
        cluster_num = cluster_data.num
    cluster_ngal_survive = np.zeros(cluster_num)
    cluster_ngal_all = np.zeros(cluster_num)
    cluster_m_i = np.zeros(cluster_num)
    for i in range(0, cluster_num):
        mass_index = cluster_data.get_cluster_mass_bin(i, mass_bins)
        cluster_m_i[i] = mass_index
        exclude_center = 0.0
        ngal_type = 'r200'
        cluster_ngal_survive[i] = cluster_data.get_ngal(i, m_infall, r_disrupt, exclude_center=exclude_center, ngal_type=ngal_type)[0]
        cluster_ngal_all[i]     = cluster_data.get_ngal(i, m_infall, np.inf,    exclude_center=exclude_center, ngal_type=ngal_type)[0]
     
        # if cluster_ngal_survive[i] == 0 or cluster_ngal_all[i] == 0:
        #     print("hmm, zeros")
        #     cluster_data.get_ngal(i, m_infall, r_disrupt, verbose=True, compact_central=True)
        #     # cluster_data.plot_cluster(i)
        #     plt.show()
    survive_mean = np.zeros(len(mass_bins)-1)
    survive_err = np.zeros(len(mass_bins)-1)
    survive_std = np.zeros(len(mass_bins)-1)
    for i in range(0, len(mass_bins)-1):
        slct = cluster_m_i == i
        survival_ratios = cluster_ngal_survive[slct]/cluster_ngal_all[slct]
        survive_mean[i] = np.nanmean(survival_ratios)
        survive_std[i]  = np.nanstd(survival_ratios)
        survive_err[i]  = survive_std[i]/np.sqrt(np.sum(slct))
        # print("{:.2e}->{:.2e}\n\tnum:{}\n\survival rate:{:.3f}+/-{:.3f}".format(mass_bins[i], mass_bins[i+1], np.sum(slct), survive_mean[i], survive_err[i]))
        # print("\t {}/{}".format(np.sum(cluster_ngal_survive[slct]), np.sum(cluster_ngal_all[slct])))
        # print("\t m_infall: {:.2e} r_disrupt: {:.2e}".format(m_infall, r_disrupt))
    return survive_mean, survive_err


def get_survival_rate(param_fname, fit, mass_bins):
    param = dtk.Param(param_fname)
    cluster_loc = param.get_string('cluster_loc')
    cluster_load_num = param.get_int('cluster_load_num')
    cluster_data = load_clusters(cluster_loc)
    # getting the Ngal with and without disruption
    model_fit_fname = "figs/"+param_fname+"/calc_likelihood_bounds.py/grid_fit_param.txt"
    model_fit = load_fit_limits(model_fit_fname)
    return get_survival_rate_cluster(cluster_data,  cluster_load_num, mass_bins, model_fit, )

def plot_luminosity_dependence_survival_single(mstars, param_base, color):
    param_fnames = [param_base.replace("@val@", str(mstar)) for mstar in mstars]
    fits = load_fit_limits_set(param_fnames)
    mass_bins_list = [ ZMR('output/'+param_fname+'/zmr_lkhd_cores.param').m_bins for param_fname in param_fnames]
    survival_infos = [ get_survival_rate(param_fname, fit, mass_bins) for param_fname, fit, mass_bins in zip(param_fnames, fits, mass_bins_list)]
    colors = ['r', 'b', 'g', 'c']
    labels_dict = {-1: "2.50~L$_*$", 
              0: "1.00~L$_*$", 
              0.5: "0.63~L$_*$", 
              1: "0.40~L$_*$"}
    labels = [ labels_dict[mstar] for mstar in mstars ]
    x_centers = dtk.bins_avg(mass_bins)
    survival_rate = []
    survival_rate_err = []
    for survival_info, mass_bins, label, color in zip(survival_infos, mass_bins_list, labels, colors):
        print(label, survival_info[0])
        plt.plot(dtk.bins_avg(mass_bins), survival_info[0],  label=label, color=color)
        plt.fill_between(dtk.bins_avg(mass_bins), survival_info[0]+survival_info[1], survival_info[0]-survival_info[1], color=color, alpha=0.3)
        survival_rate.append(survival_info[0])
        survival_rate_err.append(survival_info[1])
    return x_centers, survival_rate, survival_rate_err

def finalize_luminosity_dependence_survival_plot():
    plt.xscale('log')
    plt.xlabel(r'M$_{200m}$ [ h$^{-1}$ M$_\odot$]')
    plt.ylabel('survival fraction')
    plt.ylim([0,1])
    plt.xlim([1e14, 8e15])
    plt.legend(loc='best', framealpha=0.0)
    plt.tight_layout()
    
def plot_survival_ratio(OR,  QC):
    centers_or, survival_or, survival_err_or = OR
    centers_qc, survival_qc, survival_err_qc = QC
    plt.figure()
    colors = ['b', 'r', 'g', 'c']
    labels = ['2.5L*', 'L*', '0.63L*', '0.4L*']
    for i in range(0, 2):
        plt.plot(cneters_or, (survival_qc[i]-survival_or[i+1])/survival_or[i+1], label=labels[i+1])
    plt.axhline(0, '--k')

def plot_luminosity_dependence_survival():
    
    mstars = [-1, 0, 0.5, 1]
    # Outer Rim
    plt.figure()
    centers_or, survival_or, survival_err_or = plot_luminosity_dependence_survival_single(mstars, "params/cfn/simet/mstar@val@/mean/a3_rd.param", 'r')
    finalize_luminosity_dependence_survival_plot()
    # Continuum
    mstars = [0, 0.5, 1]
    plt.figure()
    centers_qc, survival_qc, survival_err_qc = plot_luminosity_dependence_survival_single(mstars, "params/cfn/simet/mstar@val@/mean/qc_rd.param", 'r')
    finalize_luminosity_dependence_survival_plot()

    plot_survival_ratio([centers_or, survival_or, survival_err_or], 
                        [centers_qc, survival_qc, survival_err_qc])
   
if __name__ == "__main__":
    mstars = [-1, 0, 0.5, 1]
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a2_mi.param", mstars, ['mi', 'rd', 'rm', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rd.param", mstars, ['mi', 'rd', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rm.param", mstars, ['mi', 'rm', 'x2'])
    # plot_luminosity_dependence("params/cfn/simet/mstar@val@/mean/a_rd_rm.param", mstars, ['mi', 'rd', 'rm', 'x2'])
    # plot_luminosity_dependence_parameters_all()
    plot_luminosity_dependence_survival()
    # plot_luminosity_dependent_ngal_all("params/cfn/simet/mstar0/mean/a3_rd.param")
    # plot_luminosity_dependent_ngal_all("params/cfn/simet/mstar0/crit/qc_rd.param")

    dtk.save_figs("figs/"+__file__+"/", extension=".pdf")
    plt.show()
