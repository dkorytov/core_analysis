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

def plot_1d_likelihood(label, data_bins, lkhd_data, fit_data, data_bds, log = False):
    print(",,,", label)
    plt.figure()
    plt.plot(data_bins, lkhd_data, 'x-', label='likelihood: {:.3f}->{:.3f}'.format(data_bins[data_bds[0]], data_bins[data_bds[1]]))
    ylim = plt.ylim()
    plt.ylim([0, ylim[1]])
    b1,b2 = data_bds[0], data_bds[1]+1
    plt.fill_between(data_bins[b1:b2], 0, lkhd_data[b1:b2], lw=0.0, alpha=0.3)
    plt.axvline(fit_data,c='k', ls='--')
    plt.plot([],[],'k--',label='grad descn: {:.3f}'.format(fit_data))
    max_val = np.argmax(lkhd_data)
    plt.axvline(data_bins[max_val],c='r', ls='--')
    plt.plot([],[],'r--',label='max lkhd: {:.3f}'.format(data_bins[max_val]))
    plt.grid()
    plt.legend(loc='best', framealpha=0.3)
    plt.xlabel(label)
    plt.ylabel('~ Likelihood ')
    if log:
        plt.yscale('log')


def plot_2d_likelihood(labels, data_bins, lkhd_data):
    dx = (data_bins[0][1]-data_bins[0][0])/2.0
    dy = (data_bins[1][1]-data_bins[1][0])/2.0
    plt.figure()
    cs = plt.pcolor(data_bins[0]-dx, data_bins[1]-dy, lkhd_data.T, cmap='Blues')

    #cs_lines = plt.contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, c='k')
    print(data_bins[0].shape)
    print(data_bins[1].shape)
    print(lkhd_data.shape)
    cs_lines = dtk.quick_contour(data_bins[0], data_bins[1], lkhd_data, 
                                 levels = (0.68, 0.87),
                                 colors = 'k',
                                 label = False,
                                 smoothen = False,
                                 bins_edges = False)
    # ax = plt.gca()
    # ax.axvline(fit_data[0],c='k', ls='--')
    # ax.axhline(fit_data[1],c='k', ls='--')
    
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.grid()
    cb = plt.colorbar(cs)
    # cb.add_lines(cs_lines)

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
    #fit_indx = np.searchsorted(bins,fit)
    fit_indx = np.argmax(vals)
    r_indx = fit_indx 
    l_indx = fit_indx 
    current_sum = vals[fit_indx]
    # plt.figure()
    # plt.plot(bins,vals,'-x')
    # plt.plot(bins[fit_indx],vals[fit_indx],'o')
    # ax = plt.gca()
    # i = 0
    while(True):
        print(l_indx, r_indx,end='')
        if current_sum >= limit:
            break
        if r_indx == vals.size-1:
            while(current_sum >= limit):
                current_sum += vals[l_indx]
                l_indx -= 1
        if l_indx == 0:
            while(current_sum >= limit):
                current_sum += vals[r_indx]
                r_indx += 1
        if vals[r_indx+1] > vals[l_indx-1]:
            r_indx += 1
            current_sum += vals[r_indx]
            # ax.annotate(str(i), xy=(bins[r_indx],vals[r_indx]))
        else:
            l_indx -= 1
            current_sum += vals[l_indx]
        #     ax.annotate(str(i), xy=(bins[l_indx],vals[l_indx]))
        # i +=1
        # print("->" , l_indx, r_indx, '\t', current_sum)

    # print("final: ", l_indx, r_indx)
    if l_indx< 0:
        l_indx = 0
    if r_indx == vals.size:
        r_indx -=1
    return l_indx, r_indx
            

def calc_likelihood_bounds(param_file_name):
    param = dtk.Param(param_file_name)
    lgrid_param = dtk.Param("output/"+param_file_name+"/lgrid.param")
    has_rm = param.get_bool("fit_r_merger")
    has_rd = param.get_float_list("rd_bins_info")[2]>2
    
    result = np.array(lgrid_param.get_double_list("result"))
    hfile_fit = h5py.File("output/"+param_file_name+"/fit_core_params.hdf5")
    fit_mi = hfile_fit['m_infall'].value[0]
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
    lkhd_mi_rm = np.sum(lkhd, axis=1)
    #print(np.shape(lkhd_mi))

    mi_bds = get_bounds_limits(lkhd_mi, np.log10(mi_bins), np.log10(fit_mi))
    mi_bounds = np.log10(mi_bins)[mi_bds[0]] - np.log10(fit_mi), np.log10(mi_bins)[mi_bds[1]] - np.log10(fit_mi)
    print("mi: ",np.log10(fit_mi),mi_bounds)
    if has_rd:
        fit_rd = hfile_fit['r_disrupt'].value[0]
        rd_bds = get_bounds_limits(lkhd_rd, rd_bins, fit_rd)
        rd_bounds = rd_bins[rd_bds[0]], rd_bins[rd_bds[1]]
        print("rd: ",fit_rd, rd_bounds)
    if has_rm:
        fit_rm = hfile_fit['r_merger'].value[0]
        rm_bds = get_bounds_limits(lkhd_rm, np.log10(rm_bins), np.log10(fit_rm))
        rm_bounds = np.log10(rm_bins)[rm_bds[0]] - np.log10(fit_rm), np.log10(rm_bins)[rm_bds[1]] - np.log10(fit_rm)
        print("mi: ",np.log10(fit_mi),mi_bounds)

    plot_1d_likelihood("Minfall", np.log10(mi_bins), lkhd_mi, np.log10(fit_mi), mi_bds)

    
    # plt.figure()
    # plt.plot(np.log10(mi_bins), lkhd_mi, label='likelihood')
    # ylim = plt.ylim()
    # plt.ylim([0, ylim[1]])
    # plt.fill_between(np.log10(mi_bins[mi_bds[0]:mi_bds[1]]), 0, lkhd_mi[mi_bds[0]:mi_bds[1]], lw=0.0, alpha=0.3)
    # plt.axvline(np.log10(fit_mi),c='k', ls='--')
    # plt.plot([],[],'k--',label='grad descn fit value')

    # plt.grid()
    # plt.legend(loc='best')
    # plt.xlabel('M$_{infall}$ [h$^{-1}$M$_\odot$]')
    # plt.ylabel('Likelihood ')

    if has_rd:
        plot_1d_likelihood('R$_{disrupt}$ [h$^{-1}$Mpc]',rd_bins, lkhd_rd, fit_rd, rd_bds)

    if has_rm:
        plot_1d_likelihood('R$_{merger}$ [h$^{-1}$Mpc]',rm_bins, lkhd_rm, fit_rm, rm_bds)

    # plt.figure()
    # plt.plot(rm_bins, lkhd_rm)

    # plt.figure()

    # plt.grid()
    # plt.xlabel('log10(M_infall/h^-1 Msun)')
    # plt.ylabel('R_disrupt [h^-1 Mpc]')
    if has_rd:
        plot_2d_likelihood(('M_infall', 'R_merger'), (np.log10(mi_bins), rd_bins), lkhd_mi_rd)

    if has_rm:
        plot_2d_likelihood(('M_infall', 'R_merger'), (np.log10(mi_bins), rm_bins), lkhd_mi_rm)
    # plt.figure()
    # cs = plt.pcolor(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, cmap='Blues')

    # #cs_lines = plt.contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, c='k')
    # cs_lines = dtk.quick_contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd, 
    #                              levels = (0.68, 0.87),
    #                              colors = 'k',
    #                              label = False,
    #                              smoothen = True,
    #                              bins_edges = False)
                      

    # plt.xlabel('log10(M$_{infall}$ h$^{-1}$ Msun)')
    # plt.ylabel('R$_{disrupt}$ [h$^{-1} Mpc]')
    # plt.grid()
    # cb = plt.colorbar(cs)
    # cb.add_lines(cs_lines)

    # plt.figure()
    # plt.title('Confidence Contours')

    # plt.xlabel('M$_{infall}$ [h$^{-1}$M$_\odot$]')
    # plt.ylabel('R$_{disrupt}$ [h$^{-1}$Mpc]')
    # plt.grid()
    # plt.tight_layout()

    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/')
    txt_file = file("figs/"+param_file_name+"/"+__file__+"/fit_param.txt", 'a')
    txt_file.write("mi_limits \t{}\t{}\n".format(mi_bins[mi_bds[0]],mi_bins[mi_bds[1]]))
    if has_rd:
        txt_file.write("rd_limits \t{}\t{}\n".format(rd_bins[rd_bds[0]],rd_bins[rd_bds[1]]))
    if has_rm: 
        txt_file.write("rm_limits \t{}\t{}\n".format(rm_bins[rm_bds[0]],rm_bins[rm_bds[1]]))


    
    plt.show()


def write_fit_param(param_file):
    fname = "output/"+param_file+"/fit_core_params.hdf5"
    print(fname)
    grad_descent_hfile = h5py.File(fname,'r')
    fname = "figs/"+param_file+"/"+__file__+"/fit_param.txt"
    dtk.ensure_dir(fname)
    txt_file = file(fname, 'w')
    for key in grad_descent_hfile.keys():
        txt_file.write(key+"\t"+str(grad_descent_hfile[key].value[0])+"\n")
    grad_descent_hfile.close()
    txt_file.close()

if __name__ == "__main__":
    write_fit_param(sys.argv[1])
    calc_likelihood_bounds(sys.argv[1])
    
