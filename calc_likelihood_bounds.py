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
from core_fit2_util import *
from catalog_reader import Catalog,frag_to_real
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize

def renormalize(vals, bins):
    print("norm", np.shape(vals))
    width = bins[1]-bins[0]
    print(np.shape(width))
    vals2 = vals/width
    print(np.shape(vals2))
    a = vals2/np.sum(vals2)
    print(np.shape(a))
    return a


def get_bounds_limits(vals, bins, fit, limit=0.67):
    tot = np.sum(vals)
    vals = vals/tot
    fit_indx = np.searchsorted(bins,fit)
    r_indx = fit_indx +1
    l_indx = fit_indx -1
    current_sum = vals[fit_indx]
    while(True):
        if vals[r_indx] > vals[l_indx]:
            current_sum += vals[r_indx]
            r_indx += 1
        else:
            current_sum += vals[l_indx]
            l_indx -= 1
        if current_sum >= limit:
            break
        if r_indx == vals.size-1:
            while(current_sum >= limit):
                current_sum += vals[l_indx]
                l_indx -= 1
        if l_indx == vals.size-1:
            while(current_sum >= limit):
                current_sum += vals[r_indx]
                r_indx += 1
    return l_indx, r_indx
            

def calc_likelihood_bounds(param_file_name):

    lgrid_param = dtk.Param("output/"+param_file_name+"/lgrid.param")
    result = np.array(lgrid_param.get_double_list("result"))
    hfile_fit = h5py.File("output/"+param_file_name+"/fit_core_params.hdf5")
    fit_mi = hfile_fit['m_infall'].value[0]
    fit_rd = hfile_fit['r_disrupt'].value[0]
    nan_slct = np.isnan(result)
    result[nan_slct] = np.ones(np.sum(nan_slct))*1000000
    mi_bins = np.array(lgrid_param.get_float_list("mi_bins"))
    rd_bins = np.array(lgrid_param.get_float_list("rd_bins"))
    rm_bins = np.array(lgrid_param.get_float_list("rm_bins"))
    
    result2 = result.reshape((mi_bins.size,rd_bins.size,rm_bins.size))
    
    #print(np.min(result2), np.max(result2))
    lkhd = np.exp(-(result2-np.min(result2))/2.0)
    
    lkhd_mi = np.sum(lkhd, axis=(1,2))
    lkhd_rd = np.sum(lkhd, axis=(0,2))
    lkhd_rm = np.sum(lkhd, axis=(0,1))
    # print(np.shape(lkhd_mi))
    # lkhd_mi = renormalize(lkhd_mi, np.log10(mi_bins))
    # lkhd_mi = renormalize(lkhd_rd, rd_bins)
    # lkhd_mi = renormalize(lkhd_rm, rm_bins)
    lkhd_mi_rd = np.sum(lkhd, axis=2)
    #print(np.shape(lkhd_mi))
    mi_bds = get_bounds_limits(lkhd_mi, np.log10(mi_bins), np.log10(fit_mi))
    rd_bds = get_bounds_limits(lkhd_rd, rd_bins, fit_rd)
    mi_bounds = np.log10(mi_bins)[mi_bds[0]] - np.log10(fit_mi), np.log10(mi_bins)[mi_bds[1]] - np.log10(fit_mi)
    rd_bounds = rd_bins[rd_bds[0]], rd_bins[rd_bds[1]]
    print("mi: ",np.log10(fit_mi),mi_bounds)
    print("rd: ",fit_rd, rd_bounds)

    plt.figure()
    plt.plot(np.log10(mi_bins), lkhd_mi, label='likelihood')
    plt.axvline(np.log10(fit_mi),c='k', ls='--')
    ylim = plt.ylim()
    plt.ylim([0, ylim[1]])
    plt.fill_between(np.log10(mi_bins[mi_bds[0]:mi_bds[1]]), 0, lkhd_mi[mi_bds[0]:mi_bds[1]], lw=0.0, alpha=0.3)
    plt.plot([],[],'k--',label='best fit value')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('M$_{infall}$ [h$^{-1}$M$_\odot$]')
    plt.ylabel('Likelihood ')

    plt.figure()
    plt.plot(rd_bins, lkhd_rd, label='likelihood')
    ylim = plt.ylim()
    plt.ylim([0, ylim[1]])
    plt.fill_between(rd_bins[rd_bds[0]:rd_bds[1]], 0, lkhd_rd[rd_bds[0]:rd_bds[1]], lw=0.0, alpha=0.3)
    plt.axvline(fit_rd, c='k', ls ='--')
    plt.plot([],[],'k--',label='best fit value')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('R$_{disrupt}$ [h$^{-1}$Mpc]')
    plt.ylabel('Likelihood ')


    plt.figure()
    plt.plot(rm_bins, lkhd_rm)

    plt.figure()
    plt.contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, c='k')
    plt.grid()
    plt.xlabel('log10(M_infall/h^-1 Msun)')
    plt.ylabel('R_disrupt [h^-1 Mpc]')

    plt.figure()
    plt.pcolor(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, cmap='Blues')
    plt.xlabel('log10(M$_infall$/h$^{-1}$ Msun)')
    plt.ylabel('R_disrupt [h$^{-1} Mpc]')
    plt.grid()
    plt.colorbar()

    plt.figure()
    plt.title('Confidence Contours')
    dtk.quick_contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd, 
                      levels = (0.68, 0.87),
                      colors = 'k',
                      label = False,
                      bins_edges = False)

    plt.xlabel('M$_{infall}$ [h$^{-1}$M$_\odot$]')
    plt.ylabel('R$_{disrupt}$ [h$^{-1}$Mpc]')
    plt.xlim([11.9, 12.3])
    plt.ylim([0.02, 0.07])
    plt.grid()

    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/')
    plt.show()

    
    
if __name__ == "__main__":
    calc_likelihood_bounds(sys.argv[1])
