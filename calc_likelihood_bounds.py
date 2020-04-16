#!/usr/bin/env python3

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

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)


def plot_1d_likelihood(label, data_bins, lkhd_data, fit_data, data_bds, log = False):
    data_bins = 10**data_bins
    plt.figure()
    plt.plot(data_bins, lkhd_data, '-')#, label='likelihood: {:.3f}->{:.3f}'.format(data_bins[data_bds[0]], data_bins[data_bds[1]]))
    ylim = plt.ylim()
    plt.ylim([0, ylim[1]])
    b1,b2 = data_bds[0], data_bds[1]+1
    plt.fill_between(data_bins[b1:b2], 0, lkhd_data[b1:b2], lw=0.0, alpha=0.3)
    # plt.plot([],[],'k--',label='grad descn: {:.3f}'.format(fit_data))
    max_val = np.argmax(lkhd_data)
    plt.axvline(data_bins[max_val],c='r', ls='--')
    # plt.plot([],[],'r--',label='max lkhd: {:.3f}'.format(data_bins[max_val]))
    central_value = np.log10(data_bins[max_val])
    upper_lim = np.log10(data_bins[data_bds[1]])
    lower_lim = np.log10(data_bins[data_bds[0]])
    title_str ="{} = ${:.3f}^{{+{:.3f}}}_{{-{:.3f}}}$".format(label, central_value, upper_lim-central_value, central_value-lower_lim,)
    print(title_str)
    plt.title(title_str)
    plt.grid()
    plt.legend(loc='best', framealpha=0.3)
    plt.xlabel(label)
    plt.ylabel('~ Likelihood ')
    # plt.xscale('log')
    
    if log:
        plt.yscale('log')


def plot_2d_likelihood(labels, data_bins, lkhd_data):
    dx = (data_bins[0][1]-data_bins[0][0])/2.0
    dy = (data_bins[1][1]-data_bins[1][0])/2.0
    plt.figure()
    cs = plt.pcolor(data_bins[0]-dx, data_bins[1]-dy, lkhd_data.T, cmap='Blues')

    #cs_lines = plt.contour(np.log10(mi_bins), rd_bins, lkhd_mi_rd.T, c='k')
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
    vals2 = vals/width
    a = vals2/np.sum(vals2)
    return a


def interpolate_fine_grain(vals, bins, num, fit_index = None):
    new_bins = np.linspace(np.min(bins), np.max(bins), num)
    new_vals = np.interp(new_bins, bins, vals)
    if fit_index is None:
        return new_vals, new_bins, None
    else:
        new_fit_index = np.searchsorted(new_bins, bins[fit_index])
        return new_vals, new_bins, new_fit_index


def get_bounds_limits(vals, bins, fit, limit=0.67, fine_grain=None, fit_index=None):
    if fine_grain is not None:
        vals_new, bins, fit_index = interpolate_fine_grain(vals, bins, fine_grain, fit_index=fit_index)
        vals = vals_new
    tot = np.sum(vals)
    vals = vals/tot
    fit_indx = np.searchsorted(bins,fit)
    if fit_index is None:
        fit_indx = np.argmax(vals)
    else:
        fit_indx= fit_index
        print("fit index: ", fit_indx, np.argmax(vals))
    r_indx = fit_indx 
    l_indx = fit_indx 
    current_sum = vals[fit_indx]
    # plt.figure()
    # plt.plot(bins,vals,'-x')
    # plt.plot(bins[fit_indx],vals[fit_indx],'o')
    # ax = plt.gca()
    # i = 0
    while(True):
        # print(l_indx, r_indx, current_sum)
        if current_sum >= limit:
            break
        if r_indx == vals.size-1:
            # print("r_indx_max")
            # print(current_sum,"?>=", limit)
            while(current_sum <= limit):
                current_sum += vals[l_indx]
                l_indx -= 1
                # print(current_sum,"?>=", limit)
            break
        if l_indx == 0:
            while(current_sum <= limit):
                current_sum += vals[r_indx]
                r_indx += 1
            break
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

    # print("pre-final: ", l_indx, r_indx)
    if l_indx< 0:
        l_indx = 0
    if r_indx == vals.size:
        r_indx -=1
    # print("final: ", l_indx, r_indx)
    # plt.show()
    # exit()
    if fine_grain is None:
        return l_indx, r_indx
    else:
        return l_indx, r_indx, vals_new, bins
            

def get_1d_axis_sum(dim, dims):
    result = []
    for i in dims:
        pass


def corner_plot(labels, grid_dic = None, mcmc_dic = None,
                expected_comov_abundance = None, core_loc=None):
    """Plot a corner plot for either from the likelihood calculated on a
    grid and/or from an MCMC.

    """
    size = len(labels)
    #If there is anything to plot at all
    if grid_dic is None and mcmc_dic is None and cost_dic is None:
        return 
    #Make sure the plotting options are compatible 
    if grid_dic is not None and mcmc_dic is not None:
        assert grid_dic['cost'] == None, "cost fucntion plotting cannot work with mcmc plotting"
     
    # Create the figure
    fig, axs = plt.subplots(size, size,  sharex = 'col', figsize=(10,8), squeeze=True)
    #Plotting
    if grid_dic is not None and mcmc_dic is not None:
        corner_plot_mcmc(labels, mcmc_dic['mcmc_loc'], fig, axs, colors = ['r', 'r', 'r'], plot_hist = False, alpha = 0.1)
        corner_plot_grid(labels, grid_dic['bins'], grid_dic['lkhd'], fig, axs, cost=grid_dic['cost'], colors = ['b', 'b', 'b'], plot_hist=True, alpha=0.1)

    elif grid_dic is not None:
        corner_plot_grid(labels, grid_dic['bins'], grid_dic['lkhd'], fig, axs, cost=grid_dic['cost'])
    elif mcmc_dic is not None:
        corner_plot_mcmc(labels, mcmc_dic['mcmc_loc'], fig, axs)
    #Adding axis labels
    for i in range(0,size):
        ax = axs[size-1,i]
        ax.set_xlabel(labels[i])
    for i in range(1, size):
        ax = axs[i,0]
        ax.set_ylabel(labels[i])
    #Formatting
    fig.tight_layout()
    fig.subplots_adjust(hspace=0,wspace=0)
    plt.ticklabel_format(useOffset=False, style='plain')
    # Get rid of upper right corner of subplots
    for i in range(size):
        for j in range(size):
            if i < j: 
                axs[i][j].set_visible(False)
            if i == j and grid_dic:
                if grid_dic['cost'] is not None:
                    pass
                    # ylim = axs[i][j].get_ylim()
                    # axs[i][j].set_ylim([50,1e4])
                # ylim=[0.8*ylim[0], ylim[1]]
                # axs[i][j].set_ylim(ylim)
                # yticks = axs[i][j].get_yticklabels()
                # if len(yticks) > 0:
                #     # plt.setp(axs[i][j].get_yticklabels()[-1], visible=False)
                #     plt.setp(yticks[-1], visible=False)
                    # print(help(yticks[-1]))
                    # yticks[-1] = ""
                    # axs[i][j].set_yticklabels(yticks)
                
    # Add the abundance line to rd/disrupt plot
    if expected_comov_abundance is not None:
        disrupt_index = -1
        infall_index = -1
        for i in range(0,len(labels)):
            if 'disrupt' in labels[i]:
                disrupt_index = i
            elif 'infall' in labels[i]:
                infall_index = i
        if disrupt_index != -1 and infall_index != -1:
            hfile_fname = 'tmp_hdf5/{}/abundance={}.hdf5'.format(core_loc, expected_comov_abundance)
            print(hfile_fname)
            exists = os.path.isfile(hfile_fname)
            if exists:
                hfile = h5py.File(hfile_fname, 'r')
                abund_infall_mass = hfile['abund_infall_mass'][()]
                abund_radius = hfile['abund_radius'][()]
                ax= axs[infall_index][disrupt_index]
                ax= axs[disrupt_index][infall_index]
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                ax.plot(np.log10(abund_infall_mass), abund_radius, '--k', lw=2.0)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
        

def corner_plot_mcmc(labels, mcmc_loc, fig, axs, colors =['tab:blue', 'k', 'tab:red'], plot_hist = True, alpha = 0.3):
    mcmc_m_i  = dtk.gio_read(mcmc_loc,"mcmc_mass_infall")
    mcmc_r_d  = dtk.gio_read(mcmc_loc,"mcmc_r_disrupt")
    mcmc_id   = dtk.gio_read(mcmc_loc,"mcmc_walker_id")
    mcmc_step = dtk.gio_read(mcmc_loc,"mcmc_walker_step")
    mcmc_val  = dtk.gio_read(mcmc_loc,"mcmc_value")
    size = len(labels)
    data = []
    slct = mcmc_step > np.max(mcmc_step)/2.0
    data.append(np.log10(mcmc_m_i[slct]))
    data.append(mcmc_r_d[slct])
    best_fit_indx = np.argmin(mcmc_val[slct])
    #Diagonal Covariance 
    for i in range(0, len(labels)):
        ax = axs[i][i]
        h, xbins = np.histogram(data[i], bins = 50, density=True)
        ax.plot(dtk.bins_avg(xbins), h, c=colors[0])

        # Calculate 1 simga limits
        b1, b2 = get_bounds_limits(h, xbins, 0)
        b2 += 0
        b1 += 0
        b1a = b1-1
        ylim = ax.get_ylim()
        ax.fill_between(dtk.bins_avg(bbins[b1a:b2]), 0, h[b1:b2], lw=0.0, alpha=alpha, color=colors[0])
        ax.set_xlim(np.min(xbins), np.max(xbins))
        ax.set_ylim(ylim) #restore old ylims before fill_bewteen
        ax.axvline(data[i][best_fit_indx], c=colors[2], ls='--')
    for i in range(size):
        for j in range(size):
            if i <= j:
                continue
            ax = axs[i][j]
            h, xbins,ybins = np.histogram2d(data[j], data[i], bins = 50)
            if plot_hist:
                ax.pcolor(xbins,ybins, h.T, cmap='Greys')
            dtk.quick_contour(xbins, ybins, h, ax=ax,
                              levels = (0.68, 0.87),
                              colors = colors[1],
                              label = False,
                              smoothen = False,
                              bins_edges = True)
            ax.axvline(data[j][best_fit_indx], c=colors[2], ls='--')
            ax.axhline(data[i][best_fit_indx], c=colors[2], ls='--')
            plt.sca(ax)
            plt.xticks(rotation=45)

    
def corner_plot_grid(labels, bins, lkhd, fig, axs, cost=None, colors=['tab:blue', 'k', 'tab:red'], plot_hist = True, alpha = 0.3):
    assert len(labels) == len(bins), "label count [{}] doesn't match bin count[{}]".format(len(labels), len(bins))
    assert len(labels) == len(np.shape(lkhd)), "Likelihood dims [{}] don't match label count [{}]".format(len(np.shape(lkhd)), len(labels))
    for i in range(0, len(labels)):
        assert len(bins[i]) == np.shape(lkhd)[i], "Likelihood matrix length [{}] on dim [{}, {}] doesn't match bin length [{}]. Lkhd matrix: {}".format(np.shape(lkhd[i])[i], i, labels[i], len(bins[i]), np.shape(lkhd))
    size = len(labels)
    if cost is None:
        # fig.suptitle("Likelihood")
        pass
    else:
        min_cost = np.min(cost)
        # fig.suptitle("           Fit Cost = {:.1f} X^2_red={:.2f}".format(min_cost, min_cost/(5*15-len(labels))))
    axs_list = np.arange(0,size)
    limits = []
    max_lkhd = [] #Maximum likelihood value
    max_lkhd_indx = np.unravel_index(np.argmax(lkhd, axis=None), lkhd.shape)
    # HERE
    # max_lkhd_limits = [] #maxiumum likelihood 1-sigma indexs
    # max_lkhd_limits_val = []#maxiumum likelihood 1-sigma values
    for i in range(size):
        ax = axs[i][i]
        lkhd_1d = np.sum(lkhd, axis=tuple(np.delete(axs_list, i)))
        dx = bins[i][1]-bins[i][0]
        lkhd_1d = lkhd_1d/dx/np.sum(lkhd_1d)
        limits.append((np.min(bins[i]), np.max(bins[i])))
        # max_lkhd.append(bins[i][np.argmax(lkhd_1d)])
        max_lkhd.append(bins[i][max_lkhd_indx[i]])
        # max_lkhd_limits.append(get_bounds_limits(lkhd_1d, bins[i], max_lkhd[i]))
        # max_lkhd_limits_val.append((bins[i][max_lkhd_limits[i][0]], bins[i][max_lkhd_limits[i][1]]))
        b1, b2, vals, bbins = get_bounds_limits(lkhd_1d, bins[i], max_lkhd[i], fine_grain=5000, fit_index=max_lkhd_indx[i])
        # b1, b2 = max_lkhd_limits[i]
        b2 = b2 + 1
        if b2 == len(bbins):
            b2 = len(bbins)-1
    
        if cost is None:
            ax.plot(bins[i], lkhd_1d, c=colors[0])
            # ylim = ax.get_ylim()
            ax.fill_between(bbins[b1:b2], 0, vals[b1:b2], lw=0.0, alpha=alpha, color=colors[0])
        else:
            t1 = np.amin(cost, axis=tuple(np.delete(axs_list, i)))
            ax.plot(bins[i], t1.flatten())
        ax.axvline(max_lkhd[i], c=colors[2],ls='--')
        ax.set_xlim(limits[i])
        # if cost is None:
        #     ax.set_ylim(ylim)
        plt.sca(ax)
        plt.xticks(rotation=45)
        lim_plus = bbins[b2] - max_lkhd[i]
        lim_minus = max_lkhd[i] -bbins[b1]
        if "disrupt" not in labels[i]:
            if cost is None:
                test = "$\mathrm{{ {:.3f}^{{ +{:.3f} }}_{{ -{:.3f} }} }}$".format( max_lkhd[i], lim_plus, lim_minus)
            else:
                test = "$\sim$ $\mathrm{{ {:.2f} }}$".format( max_lkhd[i])
        else: #R_disrupt variable
            if cost is None:
                test = "$\mathrm{{ {:.2f}^{{ +{:.2f} }}_{{ -{:.2f} }} }}$ [h$^{{-1}}$kpc]".format(max_lkhd[i], lim_plus, lim_minus)
            else:
                test = "$\sim$ $\mathrm{{ {:.0f} }}$ [h$^{{-1}}$kpc]".format(max_lkhd[i])


        ax.set_title(test)

        if i == 0:
            if cost is None:
                ax.set_ylabel('Likelihood')
            else:
                ax.set_ylabel('-Log Likelihood')
        if cost is not None:
            ax.set_yscale('log')
        ax.yaxis.tick_right()
        ax.set_yticks([])
        ax.minorticks_off()
    # plt.show()
    for i in range(size):
        for j in range(size):
            if i == j:
                continue
            if i < j:
                continue
            # Only the right most subplots get y axis ticks & numbers
            if j != 0:
                axs[i][j].set_yticks([])

            lkhd_2d = np.sum(lkhd, axis=tuple(np.delete(axs_list, (i,j))))
            if cost is not None:
                cost_2d = np.amin(cost, axis=tuple(np.delete(axs_list, (i,j))))
            ax = axs[i][j] 
            dx = bins[j][1]-bins[j][0]
            dy = bins[i][1]-bins[i][0]
            if plot_hist:
                # add an extra column/row for xbins & ybins for pcolormesh
                # if we don't do it, it will drop the last data column from the plot
                # and the plot will have white stripes
                xbins = np.zeros(len(bins[j])+1)
                xbins[0:-1] = bins[j]
                xbins[-1] = bins[j][-1]+dx
                xbins -= dx/2.0
                
                ybins = np.zeros(len(bins[i])+1)
                ybins[0:-1] = bins[i]
                ybins[-1] = bins[i][-1]+dy
                ybins -= dy/2.0
                if cost is not None:
                    # ax.pcolormesh(bins[j]-dx/2, bins[i]-dy/2, cost_2d.T, cmap='nipy_spectral_r', norm=clr.LogNorm())
                    ax.pcolormesh(xbins, ybins, cost_2d.T, cmap='nipy_spectral_r', norm=clr.LogNorm())
                else:
                    # ax.pcolormesh(bins[j]-dx/2, bins[i]-dy/2, lkhd_2d.T, cmap='Greys')
                    ax.pcolormesh(xbins, ybins, lkhd_2d.T, cmap='Greys')

            dtk.quick_contour(bins[j], bins[i], lkhd_2d,
                              ax =ax,
                              levels=(0.68, 0.87),
                              colors = colors[1],
                              label = False,
                              smoothen = False,
                              bins_edges=False,)
            ax.axvline(max_lkhd[j], c=colors[2], ls='--')
            ax.axhline(max_lkhd[i], c=colors[2], ls='--')
            plt.sca(ax)
            plt.xticks(rotation=45)


def calc_likelihood_bounds(param_file_name):
    param = dtk.Param(param_file_name)
    expected_comov_abundance = param.get_float('expected_comov_abundance')
    core_loc = param.get_string('core_loc')
    step = param.get_int('step')
    core_loc = core_loc.replace("${step}", str(step))
    lgrid_param = dtk.Param("output/"+param_file_name+"/lgrid.param")
    has_rm = param.get_bool("fit_r_merger")
    has_rd = param.get_float_list("rd_bins_info")[2]>2
    
 
    result = np.array(lgrid_param.get_double_list("result"))
    hfile_fit = h5py.File("output/"+param_file_name+"/fit_core_params.hdf5")

    nan_slct = np.isnan(result)
    result[nan_slct] = np.ones(np.sum(nan_slct))*1000000
    mi_bins = np.array(lgrid_param.get_float_list("mi_bins"))
    rd_bins = np.array(lgrid_param.get_float_list("rd_bins"))*1000
    rm_bins = np.array(lgrid_param.get_float_list("rm_bins"))

    result2 = result.reshape((mi_bins.size,rd_bins.size,rm_bins.size))
    
    #print(np.min(result2), np.max(result2))
    lkhd = np.exp(-(result2-np.min(result2)))
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
    max_lkhd = np.unravel_index(np.argmax(lkhd, axis=None), lkhd.shape)
    fit_mi = hfile_fit['m_infall'][()][0]
    fit_mi_bds_lwr, fit_mi_bds_upr, fit_lkhd_mi, fit_mi_bins = get_bounds_limits(lkhd_mi, np.log10(mi_bins), np.log10(fit_mi),fine_grain=5000, fit_index = max_lkhd[0])
    if has_rd:
        fit_rd = hfile_fit['r_disrupt'][()][0]
        fit_rd_bds_lwr, fit_rd_bds_upr, _, fit_rd_bins = get_bounds_limits(lkhd_rd, rd_bins, fit_rd,fine_grain=5000, fit_index=max_lkhd[1])
    if has_rm:
        fit_rm = hfile_fit['r_merger'][()][0]
        fit_rm_bds_lwr, fit_rm_bds_upr, _, fit_rm_bins = get_bounds_limits(lkhd_rm, rm_bins, fit_rm,fine_grain=5000, fit_index=max_lkhd[2])
    if has_rd and not has_rm:
        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$',
                     'R$_{\mathrm{disrupt}}$ [h$^{-1}$kpc]'], grid_dic = {'bins':
                                                            [np.log10(mi_bins), rd_bins], 'lkhd': np.sum(lkhd, axis=2),
                                                            'cost': None},
                    expected_comov_abundance=expected_comov_abundance, core_loc=core_loc)
        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$', 'R$_{\mathrm{disrupt}}$ [h$^{-1}$kpc]'],
                    grid_dic = {'bins': [np.log10(mi_bins), rd_bins], 
                                'lkhd': np.sum(lkhd, axis=2),
                                'cost': np.sum(result2, axis=2)},
                    expected_comov_abundance=expected_comov_abundance, core_loc = core_loc)
        # corner_plot([r'M$_{\mathrm{infall}}$', 'R$_{\mathrm{disrupt}}$'], mcmc_dic = {'mcmc_loc': "output/{}/mcmc.gio".format(param_file_name)})
        # corner_plot([r'M$_{\mathrm{infall}}$', 'R$_{\mathrm{disrupt}}$'], grid_dic = {'bins': [np.log10(mi_bins), rd_bins], 'lkhd': np.sum(lkhd, axis=2), 'cost': None}, mcmc_dic = {'mcmc_loc': "output/{}/mcmc.gio".format(param_file_name)})
    if has_rm and not has_rd:
        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$',
                     'R$_{\mathrm{merger}}$ [h$^{-1}$Mpc]'], grid_dic = {'bins':
                                                           [np.log10(mi_bins), rm_bins], 'lkhd': np.sum(lkhd, axis=1),
                                                           'cost': None},
                    expected_comov_abundance=expected_comov_abundance, core_loc=core_loc)
        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$',
                     'R$_{\mathrm{merger}}$ [h$^{-1}$Mpc]'], grid_dic = {'bins':
                                                           [np.log10(mi_bins), rm_bins], 'lkhd': np.sum(lkhd, axis=1),
                                                           'cost': np.sum(result2,
                                                                          axis=1)},expected_comov_abundance=expected_comov_abundance, core_loc=core_loc)
    if has_rm and has_rd:

        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$', 'R$_{\mathrm{disrupt}}$ [h$^{-1}$kpc]', 'R$_{\mathrm{merge}}$ [h$^{-1}$Mpc]', ], grid_dic  = {'bins':[np.log10(mi_bins), rd_bins, rm_bins], 'lkhd': lkhd, 'cost':None}, expected_comov_abundance=expected_comov_abundance)
        corner_plot([r'log$_{10}$M$_{\mathrm{infall}}$/h$^{-1}$M$_\odot$', 'R$_{\mathrm{disrupt}}$ [h$^{-1}$kpc]', 'R$_{\mathrm{merge}}$ [h$^{-1}$Mpc]', ], grid_dic  = {'bins':[np.log10(mi_bins), rd_bins, rm_bins], 'lkhd': lkhd, 'cost':result2}, expected_comov_abundance=expected_comov_abundance)
        #corner_plot([r'M$_{\mathrm{infall}}$', 'R$_{\mathrm{disrupt}}$', 'R$_{\mathrm{merge}}$', ], [np.log10(mi_bins), rd_bins, rm_bins], lkhd, cost = result2)
    if not has_rm and not has_rd:
        # plot_1d_likelihood("M$_{infall}$", mi_bins, lkhd_mi, fit_mi, mi_bds, log=True);
        plot_1d_likelihood("M$_{infall}$  [h$^{-1}$M$_\odot$]", fit_mi_bins, fit_lkhd_mi, fit_mi, [fit_mi_bds_lwr, fit_mi_bds_upr], log=False);
        # corner_plot([r'M$_{\mathrm{infall}}$'], grid_dic = {'bins':[np.log10(mi_bins)]]
    grid_fit_fname = "figs/"+param_file_name+"/"+__file__+"/grid_fit_param.txt"
    dtk.ensure_dir(grid_fit_fname)
    txt_file = open(grid_fit_fname, 'w')
    txt_file.write("mi\t{}\n".format(np.log10(mi_bins[np.argmax(lkhd_mi)])))
    txt_file.write("mi_limits\t{}\t{}\n".format(fit_mi_bins[fit_mi_bds_lwr], fit_mi_bins[fit_mi_bds_upr]))
    param_num = 1
    if has_rd:
        txt_file.write("rd\t{}\n".format(rd_bins[np.argmax(lkhd_rd)]))
        txt_file.write("rd_limits\t{}\t{}\n".format(fit_rd_bins[fit_rd_bds_lwr], fit_rd_bins[fit_rd_bds_upr]))
        param_num += 1
    if has_rm: 
        txt_file.write("rm\t{}\n".format(rm_bins[np.argmax(lkhd_rm)]))
        # txt_file.write("rm\t{}\n".format(fit_rm))#rm_bins[np.argmax(lkhd_rm)]))
        txt_file.write("rm_limits\t{}\t{}\n".format(fit_rm_bins[fit_rm_bds_lwr], fit_rm_bins[fit_rm_bds_upr]))
        param_num += 1
    txt_file.write("cost\t{}\n".format(np.min(result2)))
    dof = (5*15)-param_num #Five mass bins w/ halos. Each halo has 15 radial bins
    txt_file.write("X_red\t{}\n".format(np.min(result2)/dof))
    

def write_fit_param(param_file):
    fname = "output/"+param_file+"/fit_core_params.hdf5"
    print(fname)
    grad_descent_hfile = h5py.File(fname,'r')
    fname = "figs/"+param_file+"/"+__file__+"/fit_param.txt"
    dtk.ensure_dir(fname)
    txt_file = file(fname, 'w')
    for key in grad_descent_hfile.keys():
        txt_file.write(key+"\t"+str(grad_descent_hfile[key][()][0])+"\n")
    grad_descent_hfile.close()
    txt_file.close()


if __name__ == "__main__":
    param_file_name = sys.argv[1]
    # write_fit_param(param_file_name)
    calc_likelihood_bounds(sys.argv[1])
    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/')
    dtk.save_figs('figs/'+param_file_name+'/'+__file__+'/', extension='.pdf')
    plt.show()    
