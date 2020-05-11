#!/usr/bin/env python3

from __future__ import print_function, division

import dtk
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np


from generate_parameter_dist import load_fit_limits

def plot_miscentering_effects(fname_pattern, distances_str, distances_val):
    param_files_rd = [fname_pattern.replace("${displacement}", distance) for distance in distances_str]
    # [print(pfd) for pfd in param_files_rd]
    rds = []
    rds_ue = []
    rds_le = []

    mis = []
    mis_ue = []
    mis_le = []
    chi2 = []
    for pfd in param_files_rd:
        fit = load_fit_limits(pfd)
        print(fit)
        mis.append(fit['mi'])
        mis_ue.append(fit['mi_limits'][0])
        mis_le.append(fit['mi_limits'][1])
        rds.append(fit['rd'])
        rds_ue.append(fit['rd_limits'][0])
        rds_le.append(fit['rd_limits'][1])
        chi2.append(fit['x2'])

    # Convert Mpc to kpc
    distances_val *= 1000
    distances_val += 1
    # plt.figure()
    # plt.plot(distances_val+0.01, mis, '-o')
    # plt.fill_between(distances_val+0.01, mis_le, mis_ue, alpha=0.3)
    # plt.xlabel('Center Displacement')
    # plt.ylabel('M$_{infall}$ [Msun/h]')
    # plt.xscale('log')
    
    # plt.figure()
    # plt.plot(distances_val+0.01, rds, '-o')
    # plt.fill_between(distances_val+0.01, rds_le, rds_ue, alpha=0.3)
    # plt.xlabel('Center Displacement')
    # plt.ylabel('R$_{disrupt}$ [kpc/h]')
    # plt.xscale('log')

    gs = gridspec.GridSpec(3, 1)
    sharex_ax = None
    fig = plt.figure(figsize=(6,10))
    ax1 = plt.subplot(gs[0, 0])
    ax1.plot(distances_val, mis, '-o')
    ax1.fill_between(distances_val, mis_le, mis_ue, alpha=0.3)
    ax1.set_xscale('log')
    ax1.set_ylabel("M$_{infall}$"+"\n"+r"[log10 h$^{-1}$ M$_\odot$]")
    ax1.xaxis.set_ticklabels([])

    # ax1.plot(distances_val+0.01, mis, '-o')
    # ax1.fill_between(distances_val+0.01, mis_le, mis_ue, alpha=0.3)

    ax2 = plt.subplot(gs[1, 0])
    ax2.plot(distances_val, rds, '-o')
    ax2.fill_between(distances_val, rds_le, rds_ue, alpha=0.3)
    ax2.set_ylabel("R$_{disrupt}$"+"\n"+r"[h$^{-1}$ kpc]")
    ax2.set_xscale('log')
    ax2.xaxis.set_ticklabels([])
    
    ax3 = plt.subplot(gs[2, 0])
    ax3.plot(distances_val, chi2, '-o')
    ax3.set_ylabel("$\widetilde\chi^2$")
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    # ax3.fill_between(distances_val+0.01, rds_le, rds_ue, alpha=0.3)
    ax3.set_xlabel('Center Displacement + 1 [h$^{-1}$ kpc]')
    gs.tight_layout(fig)


    
if __name__ == "__main__":
    AQ_pattern = "figs/params/rmba/simet/crit/mstar0/center_displace/AQ_${displacement}_rd.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    OR_pattern = "figs/params/rmba/simet/crit/mstar0/center_displace/OR_${displacement}_rd_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    OR_pattern = "figs/params/rmba/auto/make_all_OR.McClintock.high_richness.low_rez.min20.miscentering.sh/crit/mstar0/OR_rd.miscenter_${displacement}_zoom.param/calc_likelihood_bounds.py/grid_fit_param.txt"
    # Original OR values
    # distances_str = ['0.00', '0.01', '0.02', '0.022', '0.025', '0.027', '0.03','0.04', '0.05', '0.10', '0.20', '0.50', '1.00']
    # distances_val = np.array([0.00, 0.01, 0.02, 0.022, 0.025, 0.027, 0.03, 0.04, 0.05, 0.10, 0.20, 0.50, 1.00])
    # Reran OR Values
    distances_str = ['0.00', '0.01', '0.02', '0.03','0.04', '0.05', '0.075', '0.10', '0.20', '0.50', '1.00']
    distances_val = np.array([0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.20, 0.50, 1.00])

    plot_miscentering_effects(OR_pattern, distances_str, distances_val)
    dtk.save_figs('figs/'+__file__+'/', '.pdf')
    plt.show()
