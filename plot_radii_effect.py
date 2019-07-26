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

def plot_radii_effect():
    radii = [3,4,5,7,10]
    param_file_base = "figs/params/cfn/simet/mstar0/mean/rscan_r#rad#.#model#.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    param_files_mi = [param_file_base.replace("#model#", "mi").replace("#rad#", str(radi)) for radi in radii]
    param_files_rd = [param_file_base.replace("#model#", "rd").replace("#rad#", str(radi)) for radi in radii]
    print(param_files_mi)
    plot_radii_effects(param_files_mi, radii, ['mi', 'x2'])
    plot_radii_effects(param_files_rd, radii, ['mi', 'rd', 'x2'])
    plt.show()

if __name__ == "__main__":
    plot_radii_effect();
