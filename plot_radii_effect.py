#!/usr/bin/env python2.7

from __future__ import print_function, division

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np


from generate_parameter_dist import load_fit_limits

def get_param_line(param_fits, param, return_bounds=False):
    result = []
    result_lower_bound = []
    result_upper_bound = []
    print(param_fits[0].keys())
    for param_fit in param_fits:
        result.append(param_fit[param])
        if return_bounds:
            result_lower_bound.append(param_fit[param+"_lower_err"])
            result_upper_bound.append(param_fit[param+"_upper_err"])
        if param == 'rm':
            print("{} + {} - {}".format(param_fit[param], param_fit[param+"_lower_err"], param_fit[param+"_upper_err"]))
    return np.array(result), np.array(result_lower_bound), np.array(result_upper_bound)
            
def plot_radii_effects(param_files, radii, params, xlabel=None, plot_reference=None):
    parameter_labels = {'mi': r"M$_{infall}$"+"\n"+r"[log10 h$^{-1}$ M$_\odot$]",
                        'rd': r"R$_{disrupt}$"+"\n"+r"[h$^{-1}$ kpc]",
                        'rm': r"R$_{merge}$"+"\n"+r"[h$^{-1}$ Mpc]",
                        'x2': r"$\chi^2_{reduced}$",}
    param_fits = []
    for pfile in param_files:
        print(pfile)
        fit = load_fit_limits(pfile)
        param_fits.append(fit)
    size = len(params)
    gs = gridspec.GridSpec(size, 1)
    sharex_ax = None
    fig = plt.figure()
    for i, param in enumerate(params):
        if sharex_ax is None:
            ax = plt.subplot(gs[i, 0])
        else:
            ax = plt.subplot(gs[i, 0], sharex=sharex_ax)
        y, y_err_lower, y_err_upper = get_param_line(param_fits, param, param!='x2')
        ax.plot(radii, y)
        if param != 'x2':
            print(param)
            ax.fill_between(radii, y-y_err_lower, y+y_err_upper, alpha=0.3)
            if plot_reference is not None:
                ax.axhline(y[plot_reference], ls='--', color='k')
                ax.axhline(y[plot_reference]-y_err_lower[plot_reference], ls=':', color='k')
                ax.axhline(y[plot_reference]+y_err_upper[plot_reference], ls=':', color='k')
        ax.set_ylabel(parameter_labels[param])
        if param == 'x2' :
            if xlabel is not None: #x2 is always last
                ax.set_xlabel(xlabel)
        else:
            ax.xaxis.set_ticklabels([])
    gs.tight_layout(fig)
    
def plot_cluster_radial_volume():
    radii = [3,4,5,7,10]
    param_file_base = "figs/params/cfn/simet/mstar0/mean/rscan_r#rad#.#model#.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    param_files_mi = [param_file_base.replace("#model#", "mi").replace("#rad#", str(radi)) for radi in radii]
    param_files_rd = [param_file_base.replace("#model#", "rd").replace("#rad#", str(radi)) for radi in radii]
    print(param_files_mi)
    plot_radii_effects(param_files_mi, radii, ['mi', 'x2'], xlabel=r'Radius [h$^{-1}$ Mpc]')
    plot_radii_effects(param_files_rd, radii, ['mi', 'rd', 'x2'], xlabel=r'Radius [h$^{-1}$ Mpc]')

def plot_exclude_central_region():
    excludes = [0, 1, 2, 3, 4, 5]
    labels = np.array(excludes, dtype=np.float)/15.0
    param_file_base = "figs/params/cfn/simet/mstar0/mean/a_#model#_cen#exclude#.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    # param_files_mi = [param_file_base.replace("#model#", "mi").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rd = [param_file_base.replace("#model#", "rd").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rm = [param_file_base.replace("#model#", "rm").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rd_rm = [param_file_base.replace("#model#", "rd_rm").replace("#exclude#", str(exclude)) for exclude in excludes]
    # plot_radii_effects(param_files_mi, labels, ['mi', 'x2'])
    plot_radii_effects(param_files_rd, labels, ['mi', 'rd', 'x2'], xlabel='Central Exclusion Radius/R$_{200m}$', plot_reference=0)
    plot_radii_effects(param_files_rm, labels, ['mi', 'rm', 'x2'], xlabel='Central Exclusion Radius/R$_{200m}$', plot_reference=0)
    plot_radii_effects(param_files_rd_rm, labels, ['mi', 'rd', 'rm', 'x2'], xlabel='Central Exclusion Radius/R$_{200m}$', plot_reference=0)


if __name__ == "__main__":
    plot_cluster_radial_volume();
    plot_exclude_central_region()
    dtk.save_figs("figs/"+__file__+"/", extension='.pdf')
    plt.show()
