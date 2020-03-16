#!/usr/bin/env python3

from __future__ import print_function 


import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy.random as rnd
from scipy import stats
import sys
import h5py
import dtk 
import time
from calc_ngal import ClusterData

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'], })
rc('font', size=18)


def plot_radial(clusters):
    plt.figure()
    h, xbins, ybins = np.histogram2d(clusters.core_dr_r200, clusters.core_r*1000, bins=(np.linspace(0,1.5, 100), np.logspace(0, 3, 100)))
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    median = dtk.binned_median(clusters.core_dr_r200, clusters.core_r*1000.0, xbins)
    plt.plot(dtk.bins_avg(xbins), median, 'r', label='Median')
    plt.legend(loc='best', framealpha=0.0)
    plt.yscale('log')
    plt.ylabel('Core radius [h$^{-1}$ kpc]') 
    plt.xlabel('r/R$_{200c}$')

    plt.figure()
    h, xbins, ybins = np.histogram2d(clusters.core_dr_r200, clusters.core_r*1000, bins=(np.logspace(-3, 0, 100), np.logspace(0, 3, 100)))
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    median = dtk.binned_median(clusters.core_dr_r200, clusters.core_r*1000.0, xbins)
    plt.plot(dtk.bins_avg(xbins), median, 'r', label='Median')
    plt.legend(loc='best', framealpha=0.0)
    plt.ylabel('Core radius [h$^{-1}$ kpc]') 
    plt.xlabel('r/R$_{200c}$')
    plt.yscale('log')
    plt.xscale('log')

    plt.figure()
    h, xbins, ybins = np.histogram2d(clusters.core_dr_2d_r200, clusters.core_r*1000, bins=(np.linspace(0,1.5, 100), np.logspace(0, 3, 100)))
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    median = dtk.binned_median(clusters.core_dr_2d_r200, clusters.core_r*1000.0, xbins)
    plt.plot(dtk.bins_avg(xbins), median, 'r', label='Median')
    plt.legend(loc='best', framealpha=0.0)
    plt.yscale('log')
    plt.ylabel('Core radius [h$^{-1}$ kpc]') 
    plt.xlabel('r/R$_{200c}$')

    plt.figure()
    h, xbins, ybins = np.histogram2d(clusters.core_dr_2d_r200, clusters.core_r*1000, bins=(np.logspace(-3, 0, 100), np.logspace(0, 3, 100)))
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    median = dtk.binned_median(clusters.core_dr_2d_r200, clusters.core_r*1000.0, xbins)
    plt.plot(dtk.bins_avg(xbins), median, 'r', label='Median')
    plt.legend(loc='best', framealpha=0.0)
    plt.ylabel('Core radius [h$^{-1}$ kpc]') 
    plt.xlabel('r/R$_{200c}$')
    plt.yscale('log')
    plt.xscale('log')


def plot_profile(clusters):
    m_bins = np.logspace(14, 16, 9)
    cluster_mass_bin_count, _ = np.histogram(clusters.mass, bins=m_bins)
    print(cluster_mass_bin_count)
    slct_mass = (clusters.core_host_mass > m_bins[0]) & ( clusters.core_host_mass < m_bins[1])
    print(clusters.core_host_mass)
    print(np.sum(slct_mass))
    print(np.log10(m_bins))


    slct = slct_mass
    r_bins = np.linspace(0, 1, 32)
    r_bins_area = np.diff(r_bins**2*np.pi)
    r_bins_volume  = np.diff(r_bins**3*4.0/3.0*np.pi)
    
    h2d, _ = np.histogram(clusters.core_dr_2d_r200[slct], bins=r_bins)
    h3d, _ = np.histogram(clusters.core_dr_r200[slct], bins=r_bins)
    # plt.figure()
    # plt.plot(dtk.bins_avg(r_bins), h2d/r_bins_area, '-')
    # plt.yscale('log')

    # plt.figure()
    # plt.plot(dtk.bins_avg(r_bins), h3d/r_bins_volume, '-')
    # plt.yscale('log')

    num_bins = 10
    color_map = plt.cm.coolwarm
    color_map_r = plt.cm.coolwarm_r
    colors = color_map(np.linspace(0, 1, num_bins))
    sm = plt.cm.ScalarMappable(cmap=color_map, norm=clr.LogNorm(vmin=1e11, vmax=1e13))
    cluster_m_cut_count = np.sum((clusters.mass > m_bins[0]) & (clusters.mass < m_bins[1]))
    sm._A = []
    plt.figure()
    for i, m_cut in enumerate(np.logspace(11, 13, num_bins)):
        slct_m_cut = clusters.core_m > m_cut
        h, _ = np.histogram(clusters.core_dr_r200[slct & slct_m_cut], bins=r_bins)
        label = r'M$_{{infall}} >$ {:.2e}, reduction: {:.2e}'.format(m_cut, np.float(np.sum(h))/np.sum(h3d))
        plt.plot(dtk.bins_avg(r_bins), h/r_bins_volume/cluster_m_cut_count, '-', label=label, color=colors[i])
        print(label)
    cb = plt.colorbar(sm)
    cb.set_label('M$_{infall}$')
    plt.yscale('log')
    plt.ylabel('Galaxy Density [R$_{200c}$$^{-3}$]')
    plt.xlabel('r/R$_{200c}$')
    plt.tight_layout()
    print('===')
    plt.figure()
    for i, m_cut in enumerate(np.logspace(11, 13, num_bins)):
        slct_m_cut = clusters.core_m > m_cut
        h, _ = np.histogram(clusters.core_dr_2d_r200[slct & slct_m_cut], bins=r_bins)
        label = r'M$_{{infall}} >$ {:.2e}, reduction: {:.2e}'.format(m_cut, np.float(np.sum(h))/np.sum(h2d))
        plt.plot(dtk.bins_avg(r_bins), h/r_bins_area/cluster_m_cut_count, '-', label=label, color=colors[i])
        print(label)
    cb = plt.colorbar(sm)
    cb.set_label('M$_{infall}$')
    plt.yscale('log')
    plt.ylabel('Galaxy Surface Density [R$_{200c}$$^{-2}$]')
    plt.xlabel('r/R$_{200c}$')
    plt.tight_layout()


    print("\n\n")
    slct_m_infall = clusters.core_m > -1

    colors = color_map_r(np.linspace(0,1,num_bins))
    vmin = 10**-2
    vmax = 0.3
    sm = plt.cm.ScalarMappable(cmap=color_map_r, norm=clr.LogNorm(vmin=vmin, vmax=vmax))
    sm._A = []
    plt.figure()
    for i, r_cut in enumerate(np.logspace(np.log10(vmin), np.log10(vmax), num_bins)):
        slct_r_cut = clusters.core_r < r_cut
        h, _ = np.histogram(clusters.core_dr_r200[slct & slct_r_cut & slct_m_infall], bins=r_bins)
        label = r'R$_{{disrupt}} <$ {:.2e}, reduction: {:.2e}'.format(r_cut, np.float(np.sum(h))/np.sum(h2d))
        plt.plot(dtk.bins_avg(r_bins), h/r_bins_volume/cluster_m_cut_count, '-', label=label, color=colors[i])
        print(label)
    plt.yscale('log')
    cb = plt.colorbar(sm)
    cb.set_label('R$_{disrupt}$')
    plt.ylabel('Galaxy Density [R$_{200c}$$^{-3}$]')
    plt.xlabel('r/R$_{200c}$')
    plt.tight_layout()
    print('===')
    plt.figure()
    for i, r_cut in enumerate(np.logspace(np.log10(vmin), np.log10(vmax), num_bins)):
        slct_r_cut = clusters.core_r < r_cut
        h, _ = np.histogram(clusters.core_dr_2d_r200[slct & slct_r_cut & slct_m_infall], bins=r_bins)
        label = r'R$_{{disrupt}} <$ {:.2e}, reduction: {:.2e}'.format(r_cut, np.float(np.sum(h))/np.sum(h2d))
        plt.plot(dtk.bins_avg(r_bins), h/r_bins_area/cluster_m_cut_count, '-', label=label, color=colors[i])
        print(label)
    plt.yscale('log')
    cb = plt.colorbar(sm)
    cb.set_label('R$_{disrupt}$')
    plt.ylabel('Galaxy Surface Density [R$_{200c}$$^{-2}$]')
    plt.xlabel('r/R$_{200c}$')
    plt.tight_layout()

    print("\n\n")
    
    number_tot = 10
    colors = color_map(np.linspace(0,1,number_tot))
    vmin = 1e-3
    vmax = 1e-1
    sm = plt.cm.ScalarMappable(cmap=color_map, norm=clr.LogNorm(vmin=vmin, vmax=vmax))
    sm._A = []
    # f, ax = plt.subplots()
    f_2d, ax_2d = plt.subplots()
    f_3d, ax_3d = plt.subplots()
    for i, r_merger in enumerate(np.logspace(np.log10(vmin), np.log10(vmax), number_tot)):
        # effectively no mass cut and radius cut n2 merger of cores
        clusters2 = clusters.create_compute_r_merger(0, 1e12, r_merger, clusters_number = -1) 
        clusters2.set_host_halo_mass()
        clusters2.set_core_radial_distance(write_cache=False, force=True)
        slct_mass = (clusters2.core_host_mass > m_bins[0]) & ( clusters2.core_host_mass < m_bins[1])
        h_2d, _ = np.histogram(clusters2.core_dr_2d_r200[slct_mass], bins=r_bins)
        h_3d, _ = np.histogram(clusters2.core_dr_r200[slct_mass], bins=r_bins)
        label = r'R$_{{merge}} <$ {:.2e}, reduction: {:.2e}'.format(r_merger, np.float(np.sum(h_2d))/np.sum(h2d))
        print(label)
        ax_2d.plot(dtk.bins_avg(r_bins), h_2d/r_bins_area/cluster_m_cut_count, '-', label=label, color=colors[i])
        ax_3d.plot(dtk.bins_avg(r_bins), h_3d/r_bins_volume/cluster_m_cut_count, '-', label=label, color=colors[i])
    for fig, ax in [(f_2d, ax_2d), (f_3d, ax_3d)]:
        ax.set_yscale('log')
        cb = fig.colorbar(sm)
        cb.set_label('R$_{merge}$')
        if ax == ax_2d:
            ax.set_ylabel('Galaxy Surface Density [R$_{200c}$$^{-2}$]')
        else:
            ax.set_ylabel('Galaxy Density [R$_{200c}$$^{-3}$]')
        ax.set_xlabel('r/R$_{200c}$')
        fig.tight_layout()


def plot_cluster_radial_core_prop(fname):
    clusters = ClusterData()
    clusters.load_file(fname, core_host_mass=True)
    clusters.set_core_radial_distance(write_cache=False)
    # plot_radial(clusters)
    plot_profile(clusters)

    

if __name__ == "__main__":
    # plot_cluster_radial_core_prop("tmp_hdf5/clusters_AQ_M200c.hdf5")
    plot_cluster_radial_core_prop("tmp_hdf5/clusters_OR_M200c.hdf5")
    dtk.save_figs("figs/"+__file__+"/OR/", extension=".pdf")
    dtk.save_figs("figs/"+__file__+"/OR/", extension=".png")
    plt.show()
