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

def get_param_line(param_fits, param):
    result = []
    result_lower_bound = []
    result_upper_bound = []
    print(param_fits)
    for param_fit in param_fits:
        result.append(param_fit[param])
    return result
        

def plot_radii_effects(param_files, radii, params):
    param_fits = []
    for pfile in param_files:
        fit = load_fit_limits(pfile)
        param_fits.append(fit)
    size = len(params)
    gs = gridspec.GridSpec(size, 1)
    plt.figure()
    for i, param in enumerate(params):
        ax = plt.subplot(gs[i,0])
        y = get_param_line(param_fits, param)
        ax.plot(radii, y)
        ax.set_ylabel(param)
    


    # plot_radii_effects(param_fits, radii, ['mi', 'x2'])

def plot_cluster_radial_volume():
    radii = [3,4,5,7,10]
    param_file_base = "figs/params/cfn/simet/mstar0/mean/rscan_r#rad#.#model#.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    param_files_mi = [param_file_base.replace("#model#", "mi").replace("#rad#", str(radi)) for radi in radii]
    param_files_rd = [param_file_base.replace("#model#", "rd").replace("#rad#", str(radi)) for radi in radii]
    print(param_files_mi)
    plot_radii_effects(param_files_mi, radii, ['mi', 'x2'])
    plot_radii_effects(param_files_rd, radii, ['mi', 'rd', 'x2'])


def plot_exclude_central_region():
    excludes = [0, 1, 2, 3, 4, 5]
    labels = np.array(excludes, dtype=np.float)/15.0
    param_file_base = "figs/params/cfn/simet/mstar0/mean/a_#model#_cen#exclude#.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    # param_files_mi = [param_file_base.replace("#model#", "mi").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rd = [param_file_base.replace("#model#", "rd").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rm = [param_file_base.replace("#model#", "rm").replace("#exclude#", str(exclude)) for exclude in excludes]
    param_files_rd_rm = [param_file_base.replace("#model#", "rd_rm").replace("#exclude#", str(exclude)) for exclude in excludes]
    # plot_radii_effects(param_files_mi, labels, ['mi', 'x2'])
    plot_radii_effects(param_files_rd, labels, ['mi', 'rd', 'x2'])
    plot_radii_effects(param_files_rm, labels, ['mi', 'rm', 'x2'])
    plot_radii_effects(param_files_rd_rm, labels, ['mi', 'rd', 'rm', 'x2'])

if __name__ == "__main__":
    plot_cluster_radial_volume();
    plot_exclude_central_region()
    plt.show()
