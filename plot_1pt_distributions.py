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
import sys
import dtk
import h5py

def get_core_host_mass(hfile):
    ## 
    core_cluster_mass = np.zeros(hfile['cores/core_m'].size)
    cluster_size = hfile['cluster/htag'].size
    cluster_mass = hfile['cluster/sod_mass'].value
    cluster_core_offset = hfile['cluster/core_offset'].value
    cluster_core_size = hfile['cluster/core_size'].value
    for i in range(0, cluster_size):
        start = cluster_core_offset[i]
        end = cluster_core_offset[i]+cluster_core_size[i]
        core_cluster_mass[start:end] = cluster_mass[i]
        print(start, end, np.log10(cluster_mass[i]))
    return core_cluster_mass

def plot_cluster_cores2(hfile, mass_cuts, bins = 100, log=False):
    print(hfile['cores'].keys())
    core_m = np.log10(hfile['cores']['core_m'].value)
    core_r = hfile['cores']['core_r'].value
    core_c = hfile['cores']['core_is_central'].value
    core_step = hfile['cores']['core_step'].value
    slct_sat = core_step != 401
    print(np.unique(core_step))
    xbins = 10**np.linspace(np.min(np.log10(core_r)), np.max(np.log10(core_r)), 100)
    plt.figure()
    h, xbins = np.histogram(core_r[slct_sat], bins = xbins, density=False)
    plt.plot(dtk.bins_avg(xbins),h,label='All')
    for mass_cut in mass_cuts:
        slct = core_m>mass_cut
        h,xbins = np.histogram(core_r[slct & slct_sat], bins=xbins,density=False)
        plt.plot(dtk.bins_avg(xbins), h, label='M_infall > {}'.format(mass_cut))
    plt.legend(loc=2, framealpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Core Radius [Mpc/h]')
    plt.ylabel('Counts')
    plt.title("Satellite Cores")
    plt.figure()
    h, xbins = np.histogram(core_r[slct_sat], bins = xbins)
    xbins_log = np.log10(xbins)
    xbins_width = (xbins_log[1:]-xbins_log[:-1])
    plt.plot(dtk.bins_avg(xbins),h/xbins_width/np.sum(h),label='All')
    for mass_cut in mass_cuts:
        slct = core_m>mass_cut
        h,xbins = np.histogram(core_r[slct & slct_sat], bins=xbins)
        plt.plot(dtk.bins_avg(xbins), h/xbins_width/np.sum(h), label='M_infall > {}'.format(mass_cut))
    plt.legend(loc=2, framealpha=0.3)
    plt.xscale('log')
    plt.xlabel('Core Radius [Mpc/h]')
    plt.ylabel('Normalized Population Density')
    plt.title("Satellite Cores")

    plt.figure()
    h,xxbins,yybins = np.histogram2d(core_m[slct_sat], np.log10(core_r[slct_sat]), bins = 100)
    plt.pcolor(xxbins,yybins, h.T, cmap='Blues', norm=clr.LogNorm())
    plt.grid()
    plt.ylabel('Core Radius [log10 Mpc/h]')
    plt.xlabel('M_infall [log10 Msun/h]')
    plt.title("Satellite Cores")

def plot_cluster_mass(hfile, data_cluster_fname):
    plt.figure()
    mass = np.log10(hfile['cluster/sod_mass'])
    h,xbins = np.histogram(mass,bins=64, normed = True)
    plt.plot(dtk.bins_avg(xbins), h, '-o',label='Simulation [{:}]'.format(mass.size))
    print("reading the data...")
    rm_cluster_hfile = h5py.File(data_cluster_fname,'r')
    print("Data has been opened")
    i = 0
    mass_redmapper = np.zeros(26110,dtype='f4')
    while i < 26110:
        mass_redmapper[i] = rm_cluster_hfile['gal_prop'+str(i)]['mass'].value
        i += 1
    h, xbins = np.histogram(np.log10(mass_redmapper),bins=xbins, normed=True)
    plt.plot(dtk.bins_avg(xbins), h, '-o', label='RedMapper [{:}]'.format(mass_redmapper.size))
    plt.yscale('log')
    plt.legend(loc='best', framealpha=0.3)
    plt.xlabel('Halo M$_{200c}$ [Msun/h]')
    plt.ylabel('Normalized Counts')
    plt.grid()

def plot_cluster_redshift(core_redshift, data_cluster_fname):
    plt.figure()
    print("reading the data...")
    rm_cluster_hfile = h5py.File(data_cluster_fname,'r')
    print("Data has been opened")
    i = 0
    redshift_redmapper = np.zeros(26110,dtype='f4')
    while i < 26110:
        redshift_redmapper[i] = rm_cluster_hfile['gal_prop'+str(i)]['z'].value
        i += 1
    h, xbins = np.histogram(redshift_redmapper,bins=np.linspace(0,0.6,100))
    plt.plot(dtk.bins_avg(xbins), h, '-', label='RedMapper [{:}]'.format(redshift_redmapper.size))
    plt.axvline(core_redshift, ls='--', c='r', label='Core Redshift')
    plt.yscale('log')
    plt.legend(loc='best', framealpha=0.3)
    plt.xlabel('redshift')
    plt.ylabel('Counts')
    plt.grid()


def plot_1pt_distributions(param_file):
    param = dtk.Param(param_file)
    cluster_hfname = param.get_string("cluster_loc")
    zmr_hfname = param.get_string("zmrh5_loc")

    core_hfile = h5py.File(cluster_hfname, 'r')
    zmr_hfile = h5py.File(zmr_hfname, 'r')
    plot_cluster_mass(core_hfile, '/home/dkorytov/phys/Ngal_sdss/data/normal_wmap7/query_results/query_results.hdf5')
    plot_cluster_redshift(0.24, '/home/dkorytov/phys/Ngal_sdss/data/normal_wmap7/query_results/query_results.hdf5')
    plot_cluster_cores2(core_hfile,[10, 10.5, 11, 11.5, 12, 12.5, 13.0, 13.5, ])
    dtk.save_figs("figs/"+param_file+"/"+__file__+"/")
    plt.show()

if __name__ == "__main__":
    plot_1pt_distributions(sys.argv[1])
