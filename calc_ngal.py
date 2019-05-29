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
from matplotlib.colors import LogNorm
import dtk
import sys
import time
import h5py
from scipy.special import erf
from matplotlib.patches import Circle


from scipy.optimize import minimize
from generate_parameter_dist import load_fit_limits

class HODModel:
    def __init__(self, log_M_min, sigma, log_M_0, log_M_1, alpha):
        self.M_min = 10**log_M_min
        self.log_M_min = log_M_min
        self.sigma = sigma
        self.M_0 = 10**log_M_0
        self.M_1 = 10**log_M_1
        self.alpha = alpha
        # self.M_1a = 10**log_M_1a
        # self.b_ga = b_ga

    def hod_cen(self, M200m):
        return 0.5*(1.0+erf((np.log10(M200m) - self.log_M_min)/self.sigma))
        
    def hod_sat(self, M200m):
        result = self.hod_cen(M200m)*(((M200m - self.M_0)/self.M_1)**self.alpha)
        result[~np.isfinite(result)] = 0
        return result

    def hod(self, M200m):
        return self.hod_cen(M200m) + self.hod_sat(M200m)
                      
class ClusterData:
    def load_file(self, filename, treat_centrals=False, step=None, core_host_mass=False):
        hfile = h5py.File(filename, 'r')
        self.num = hfile['cluster_num'].value
        self.x =    hfile['cluster/x'].value
        self.y =    hfile['cluster/y'].value
        self.z =    hfile['cluster/z'].value
        self.rad =  hfile['cluster/sod_radius'].value
        self.mass = hfile['cluster/sod_mass'].value
        self.htag = hfile['cluster/htag'].value
        self.step = hfile['cluster/step'].value
        self.core_offset = hfile['cluster/core_offset'].value
        self.core_size = hfile['cluster/core_size'].value
        self.core_htag = hfile['cores/core_htag'].value

        # self.core_infall_htag = hfile['cores/core_infall_htag'].value
        print(hfile['cores'].keys())
        self.core_x = hfile['cores/core_x'].value
        self.core_y = hfile['cores/core_y'].value
        self.core_z = hfile['cores/core_z'].value
        self.core_r = hfile['cores/core_r'].value
        self.core_m = hfile['cores/core_m'].value
        self.core_step = hfile['cores/core_step'].value
        #self.core_infall_htag= hfile['cores/infall_htag'].value
        self.core_is_central = hfile['cores/core_is_central'].value
        if treat_centrals:
            slct_central = self.core_step == step
            self.core_r[slct_central] = 0
        self.core_host_mass = np.zeros_like(self.core_m)
        if core_host_mass:
            for i in range(0, self.num):
                start = self.core_offset[i]
                stop  = self.core_offset[i]+self.core_size[i]
                self.core_host_mass[start:stop] = self.mass[i]

    def plot_cluster(self, i):
        f, (ax1, ax2) = plt.subplots(1,2,figsize=(25,10))
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        ax1.plot(self.core_x[start:stop], self.core_y[start:stop], 'go',mfc='none',mec='g')
        
        slct_central = self.core_step[start:stop] == 401
        slct_halo_central = (self.core_step[start:stop]==401) & (self.core_htag[start:stop] == self.htag[i])
        print(np.sum(slct_halo_central))
        dx = self.x[i] - self.core_x[start:stop]
        dy = self.y[i] - self.core_y[start:stop]
        dz = self.z[i] - self.core_z[start:stop]
        dr = dx*dx + dy*dy + dz*dz
        srt = np.argsort(dr)
        
        print(self.htag[i], "{:.2e}".format(self.mass[i]))
        print(self.x[i], self.y[i], self.z[i])
        cnt = 0
        for j in range(0, np.sum(slct_halo_central)):
            print(self.core_htag[start:stop][slct_halo_central][j], self.core_step[start:stop][slct_halo_central][j], "mass:{:.2e} rad:{:.4f} dr:{:.3f} ".format(self.core_m[start:stop][slct_halo_central][j], self.core_r[start:stop][slct_halo_central][j], dr[slct_halo_central][j]))
            # print("\t", self.core_step[start:stop][srt][j]==401)
            # print("\t", self.core_htag[start:stop][srt][j]==self.htag[i])
                  #, self.core_x[start:stop][j], self.core_y[start:stop][j], self.core_z[start:stop][j])
        print("\n")

        ax1.plot(self.core_x[start:stop][slct_central],  self.core_y[start:stop][slct_central], '+r',mfc='none',mec='r', mew=2)
        ax1.plot(self.core_x[start:stop][slct_halo_central],  self.core_y[start:stop][slct_halo_central], 'x',mfc='none',mec='b', mew=2, alpha=1.0)
        slct = self.core_htag[start:stop] == self.htag[i]
        ax1.plot(self.core_x[start:stop][slct],  self.core_y[start:stop][slct], '.',mfc='g',mec='none', mew=1, alpha=0.3)
        circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax1.plot(self.x[i], self.y[i], 'xk', mew=2)
        ax1.add_artist(circle)
        ax1.grid()
        ax1.set_title("cluster[{}]".format(i))
        ax2.plot(self.core_z[start:stop], self.core_y[start:stop], 'go',mfc='none',mec='g')
        ax2.plot(self.core_z[start:stop][slct_central],  self.core_y[start:stop][slct_central], '+r',mfc='none',mec='r', mew=2)

        ax2.plot(self.core_z[start:stop][slct_halo_central],  self.core_y[start:stop][slct_halo_central], 'x',mfc='none',mec='b', mew=2, alpha=1.0)
        ax2.plot(self.core_z[start:stop][slct],  self.core_y[start:stop][slct], '.',mfc='g',mec='none', mew=1, alpha=0.3)
        circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax2.plot(self.z[i], self.y[i], 'xk', mew=2)
        ax2.add_artist(circle)
        ax2.grid()
        ax2.set_title("cores[{}]".format(self.core_size[i]))
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        
        print(self.core_m[start:stop])
        print(self.core_r[start:stop])
        plt.figure()
        h, xbins, ybins = np.histogram2d(self.core_m[start:stop], self.core_r[start:stop], bins = (np.logspace(10,15,32), np.logspace(-3,0,32)))
        plt.pcolor(xbins, ybins, h.T+0.1, cmap='Blues', norm=clr.LogNorm())
        plt.xlabel('Core Mass')
        plt.ylabel('Core Radius')
        plt.yscale('log')
        plt.xscale('log')

    def plot_find_central(self, i):
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        slct = core_m > m_infall
        slct_compact = core_r < r_disrupt
        
    def plot_fit_cluster(self, i , m_infall,  r_disrupt, r_merger=0.0, plot_3d=False):
        # if self.mass[i] > 2e14:
        #     return
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        slct = core_m > m_infall
        slct_compact = core_r < r_disrupt
        
        dx = core_x - self.x[i]
        dy = core_y - self.y[i]
        dz = core_z - self.z[i]
        dr = np.sqrt(dx*dx + dy*dy + dz*dz)
        slct_within = dr < self.rad[i]
        print(np.sum(slct))
        print(np.sum(slct & slct_compact))
        within_num =np.sum(slct_within &slct)
        within_compact_num = np.sum(slct&slct_within&slct_compact)
        within_compact_fract = within_compact_num/within_num
        print("candidates within: {}, compact candidates: {}, fract: {:.2f}".format(within_num, within_compact_num, within_compact_fract))
        if(np.sum(slct) == 0):
            return 
        print("max core mass: {:.2e}".format(np.max(core_m)))
        # srt = np.argsort(dr)
        # print("{:.2e}".format(np.max(core_m[srt][:100])))
        f, (ax1, ax2) = plt.subplots(1,2,figsize=(25,10))
        ax1.plot(core_x[slct & ~slct_compact], core_y[slct & ~slct_compact], 'go',mfc='none',mec='g')
        ax1.plot(core_x[slct &  slct_compact], core_y[slct &  slct_compact], 'go',mfc='g',mec='g')
        circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax1.add_artist(circle)
        ax1.grid()
        
        ax2.plot(core_z[slct & ~slct_compact], core_y[slct & ~slct_compact], 'go',mfc='none',mec='g')
        ax2.plot(core_z[slct &  slct_compact], core_y[slct &  slct_compact], 'go',mfc='g',mec='g')
        circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax2.add_artist(circle)
        ax2.grid()
        
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax2.set_title("cores[{}, {}, {}]".format(self.core_size[i], np.sum(slct), np.sum(slct&slct_compact)))
        ax1.set_title("cluster[{}] \n Mass {:.2e}, R200: {:.2f}".format(i, self.mass[i], self.rad[i]))

        if plot_3d:
            f = plt.figure()
            ax = f.gca(projection = "3d")
            ax.scatter(core_x[slct &  slct_compact], core_y[slct &  slct_compact], core_z[slct &  slct_compact], edgecolors='r', facecolor = 'r', s = 3, alpha=0.3)
            ax.scatter(core_x[slct & ~slct_compact], core_y[slct & ~slct_compact], core_z[slct & ~slct_compact], edgecolors='g', facecolor = 'none', s=3, alpha=0.3)
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)

            x = self.rad[i] * np.outer(np.cos(u), np.sin(v)) + self.x[i]
            y = self.rad[i] * np.outer(np.sin(u), np.sin(v)) + self.y[i]
            z = self.rad[i] * np.outer(np.ones(np.size(u)), np.cos(v)) + self.z[i]

            ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.2)
            # Create cubic bounding box to simulate equal aspect ratio
            max_range = np.array([core_x.max()-core_x.min(), core_y.max()-core_y.min(), core_z.max()-core_z.min()]).max()
            Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(core_x.max()+core_x.min())
            Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(core_y.max()+core_y.min())
            Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(core_z.max()+core_z.min())
            # Comment or uncomment following both lines to test the fake bounding box:
            for xb, yb, zb in zip(Xb, Yb, Zb):
                ax.plot([xb], [yb], [zb], 'w')

            ax.set_aspect('equal')

        # f, (ax1, ax2) = plt.subplots(1, 2, figsize = (25, 10), sharey=True)
        # x = core_x[slct &  slct_compact]
        # y = core_y[slct &  slct_compact]
        # z = core_z[slct &  slct_compact]
        # gal_x, gal_y, gal_z, gal_w = n2merger.n2merger3d(x, y, z, r_merger)
        # ax1.scatter(x, y, facecolor='g', edgecolor='g', alpha=0.3)
        # ax1.scatter(gal_x, gal_y, facecolor='r', edgecolor='r', alpha=0.3)

        # ax2.scatter(z, y, facecolor='g', edgecolor='g', alpha=0.3)
        # ax2.scatter(gal_z, gal_y, facecolor='r', edgecolor='r', alpha=0.3)

        # circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        # ax1.add_artist(circle)
        # circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        # ax2.add_artist(circle)
        
        # ax1.set_aspect('equal')
        # ax2.set_aspect('equal')
        return 

    def get_ngal(self, i, mass_cut, radius_cut, compact_central=False):
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        core_step = self.core_step[start:stop]
        slct = (core_r<radius_cut) & (core_m>mass_cut)
        if compact_central:
            slct_central = core_step == self.step
            slct[slct_central] = True
        dx = core_x[slct]-self.x[i]
        dy = core_y[slct]-self.y[i]
        dz = core_z[slct]-self.z[i]
        dr2 = dx*dx + dy*dy + dz*dz
        dr = np.sqrt(dr2)
        slct_ngal = dr2<self.rad[i]**2
        ngal = np.sum(slct_ngal)
        # plt.figure()
        # plt.plot(core_x,  core_y, 'og',mfc='none',mec='g', mew=1, alpha=0.3)
        # plt.plot(core_x[slct],  core_y[slct], 'og',mfc='g',mec='g', mew=2)
        # plt.plot(core_x[slct][slct_ngal], core_y[slct][slct_ngal], 'bx',mfc='b',mec='b', mew=2)
        # circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        # plt.gca().add_artist(circle)
        # plt.gca().set_aspect('equal')
        # plt.show()
        # See notes in src/main.cpp or test/radial_binning_test.py
        dr_out = dr[~slct_ngal]
        start_angle = 0
        end_angle = np.arcsin(self.rad[i]/dr_out)
        weight = -np.cos(end_angle) + np.cos(start_angle)
        sum_weight = np.sum(weight)
        return ngal, ngal+sum_weight

    def get_inside_core(self, mass_cut, compact_centrals=False, r_limit = None):
        core_rs = []
        core_ms = []
        core_hms = []
        core_step = []
        for i in range(0, self.num):
            start = self.core_offset[i]
            stop = start + self.core_size[i]
            core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
            core_r = self.core_r[start:stop]
            core_m = self.core_m[start:stop]
            core_step = self.core_step[start:stop]
            slct = (core_r<radius_cut) & (core_m>mass_cut)
            if compact_central:
                slct_central = core_step == self.step
                slct[slct_central] = True
            dx = core_x[slct]-self.x[i]
            dy = core_y[slct]-self.y[i]
            dz = core_z[slct]-self.z[i]
            dr2 = dx*dx + dy*dy + dz*dz
            dr = np.sqrt(dr2)
            if r_limit is None:
                slct_ngal = dr2 < self.rad[i]**2
            else:
                slct_ngal = dr2 < r_limit**2
            core_rs.append(core_r[slct_ngal])
            core_ms.append(core_m[slct_ngal])
            core_steps.append(core_step[slct_ngal])
        core_m = np.concatenate(core_ms)
        core_r = np.concatenate(core_rs)
        return core_m, core_r

class SODData:
    def load_sod(self, sod_loc, sod_hdf5):
        if sod_hdf5:
            hfile = h5py.File(sod_loc, 'r')
            self.rad = hfile['sod_halo_radius_r200m'].value
            self.mass = hfile['sod_halo_mass_m200m'].value
            self.rad_c = hfile['sod_halo_radius_r200m'].value
            self.mass_c = hfile['sod_halo_radius_r200c'].value
        else:
            self.rad = dtk.gio_read(sod_loc, 'sod_halo_radius')
            self.mass = dtk.gio_read(sod_loc, 'sod_halo_mass')
        

def get_zeng07_hod(Mr=None):
    if Mr is None:
        raise KeyError
    if Mr == -18:
        return HODModel(11.35, 0.25, 11.20, 12.40, 0.83)
    elif Mr == -18.5:
        return HODModel(11.46, 0.24, 10.59, 12.68, 0.97)
    elif Mr == -19:
        return HODModel(11.60, 0.26, 11.49, 12.83, 1.02)
    elif Mr == -19.5:
        return HODModel(11.75, 0.28, 11.69, 13.01, 1.06)
    elif Mr == -20.0:
        return HODModel(12.02, 0.26, 11.38, 13.31, 1.06)
    elif Mr == -20.5:
        return HODModel(12.30, 0.21, 11.84, 13.58, 1.12)
    elif Mr == -21.0:
        return HODModel(12.79, 0.39, 11.92, 13.94, 1.15)
    elif Mr == -21.5:
        return HODModel(13.38, 0.51, 13.94, 13.91, 1.04)
    elif Mr == -22.0:
        return HODModel(14.22, 0.77, 14.00, 14.69, 0.87)
    else:
        raise KeyError

def get_zeng07_Mrs():
    return [-18.0, -18.5, -19.0, -19.5, -20.0, -20.5, -21.0, -21.5, -22.0]

def get_fit_limits_fname(param_fname):
    result = "figs/{}/calc_likelihood_bounds.py/grid_fit_param.txt".format(param_fname)
    return result

def get_mr_limits(param_fname):
    if "/mstar0/" in param_fname:
        return -21.5, -21.0
    elif "/mstar0.5/" in param_fname:
        return -21.0, -20.5
    elif "/mstar1/" in param_fname:
        return -20.5, -20.0
    elif "/mstar-1/" in param_fname:
        return -22.5, -22.

def get_clusters(param_fname, core_host_mass = False):
    param = dtk.Param(param_fname)
    cluster_loc = param.get_string('cluster_loc')
    central = param.get_bool('force_central_as_galaxy')
    if central:
        print("It's compact central")
    clusters = ClusterData();
    clusters.load_file(cluster_loc, treat_centrals = central, step=401, core_host_mass = core_host_mass)
    return clusters, central

def get_fit_param(param_fname):
    fit_fname = get_fit_limits_fname(param_fname)
    fit_param = load_fit_limits(fit_fname)
    fit_mi = 10**fit_param['mi']
    fit_rd = fit_param['rd']*1e-3
    return fit_mi, fit_rd

def calc_ngal(param_fname):
    clusters, central = get_clusters(param_fname)
    fit_mi, fit_rd = get_fit_param(param_fname)
    # print(fit_fname)
    ngals = []
    ngals_proj = []
    masses = []
    print(clusters.num)
    for i in range(0,clusters.num):
        ngal, ngal_proj = clusters.get_ngal(i, fit_mi, fit_rd, compact_central = central)
        masses += [clusters.mass[i]]
        ngals += [ngal]
        ngals_proj+= [ngal_proj]
    masses = np.array(masses)
    ngals = np.array(ngals)
    ngals_proj = np.array(ngals_proj)
    xbins = np.logspace(14,15.5, 32)
    xbins_cen = dtk.bins_avg(xbins)
    plt.figure()
    # plt.plot(masses, ngals, ',b', alpha=0.01, label = 'Ngal')
    # plt.plot(masses, ngals_proj, ',g', alpha=0.01, label='Projected Ngal')
    ngal_avg = dtk.binned_average(masses, ngals, xbins)
    ngal_proj_avg = dtk.binned_average(masses, ngals_proj, xbins)

    plt.plot(xbins_cen, ngal_avg, 'b', lw=2, label = "Core Ngal")
    plt.plot(xbins_cen, ngal_proj_avg, 'g', lw=2, label = "Core Ngal Avg Projected")
    plt.xscale('log')
    plt.yscale('log')

    ylim = plt.ylim()
    xlim = plt.xlim()
    mr_upper, mr_lower = get_mr_limits(param_fname)
    mass = np.logspace(11,15.5,100)
    model_lower = get_zeng07_hod(mr_lower)
    hod_tot_lower = model_lower.hod(mass)
    plt.plot(mass, hod_tot_lower, '--r', lw = 1, label="Mr<{:.1f} HOD".format(mr_lower))
    model_upper = get_zeng07_hod(mr_upper)
    hod_tot_upper = model_upper.hod(mass)
    plt.plot(mass, hod_tot_upper, 'r--', lw = 1, label="Mr<{:.1f} HOD".format(mr_upper))    
    plt.plot(mass, calc_log_mean(hod_tot_lower, hod_tot_upper), 'r', label='Mr<{:.1f} HOD'.format((mr_upper+mr_lower)/2.0))
    plt.legend(loc='best')

    plt.xlim(xlim)
    plt.ylim(ylim)

    
    plt.show()

def calc_zeng07():
    mass = np.logspace(11,15.5,100)
    Mrs = get_zeng07_Mrs()
    plt.figure()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#ff7f0e', '#8c564b']
    for i, Mr in enumerate(Mrs):
        hod_model = get_zeng07_hod(Mr)
        hod_cen = hod_model.hod_cen(mass)
        hod_sat = hod_model.hod_sat(mass)
        hod_tot = hod_model.hod(mass)
        plt.plot(mass, hod_tot, lw=2, color=colors[i])
        plt.plot(mass, hod_cen, lw=1, color=colors[i], alpha=0.3, ls='--')
        plt.plot(mass, hod_sat, lw=1, color=colors[i], alpha=0.3, ls=':')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

def calc_log_mean(data1, data2):
    a = np.sqrt(data1*data2)
    return a

def calc_disruption(param_fname):
    clusters, central = get_clusters(param_fname)
    fit_mi, fit_rd = get_fit_param(param_fname)
    survival_rates = []
    t0 = time.time()
    cluster_limit = clusters.num[0]
    cluster_limit = 50000
    for i in range(0, cluster_limit):
        ngal_surv, _ = clusters.get_ngal(i, fit_mi, fit_rd, compact_central = central)
        ngal_init, _ = clusters.get_ngal(i, fit_mi, 1e10, compact_central = central)
        survival_rates += [ngal_surv/ngal_init]
        if i % 1000 == 1:
            done_amount = i/cluster_limit
            t1 = time.time()
            time_left = (t1-t0)/done_amount *(1-done_amount)
            print("{:.3f} ETA: {:.1f}".format(done_amount, time_left))
    print("total survival rate: ", np.average(survival_rates))

def wetzel09_disruption_time(Cdyn, M_infall, M_host, z):
    a = 1.0/(1.0+z)
    t_hubble = get_hubble_time(z)
    print(t_hubble)
    t_dyn = Cdyn * (M_host/M_infall) / np.log(1.0 + M_host/M_infall) * t_hubble
    return t_dyn
    
def calc_wetzel09_disruption(param_fname):
    clusters, central = get_clusters(param_fname, core_host_mass= True)
    fit_mi, fit_rd = get_fit_param(param_fname)
    stepz = dtk.StepZ(sim_name = 'AlphaQ')
    z = stepz.get_z(401)
    print(z)
    core_t_dyn = wetzel09_disruption_time(0.25, clusters.core_m, clusters.core_host_mass, z)
    print(core_t_dyn)
    plt.figure()
    lg_min = np.min(np.log10(core_t_dyn))
    lg_max = np.max(np.log10(core_t_dyn))
    print(lg_min, lg_max)
    xbins = np.logspace(lg_min, lg_max, 10)
    print(xbins)
    xbins_cen = dtk.log_bins_avg(xbins)
    h,_ = np.histogram(core_t_dyn, bins=xbins)
    plt.loglog(xbins_cen, h, )

    plt.figure()
    plt.hist2d(np.log10(clusters.core_m),
               np.log10(clusters.core_host_mass), bins=100, cmap='Blues',
               norm=clr.LogNorm())
    
    plt.figure()
    h, xbins = np.histogram(np.log10(clusters.core_host_mass/clusters.core_m), bins=100)
    xbins_cen = dtk.bins_avg(xbins)    
    plt.plot(xbins_cen, h, label='All Cores')

    plt.xlabel('log10(Mhost/Msat)')
    plt.ylabel('Counts')


    if fit_mi < 1e3:
        fit_mi = 10**fit_mi
    slct = clusters.core_m>fit_mi
    h, xbins = np.histogram(np.log10(clusters.core_host_mass[slct]/clusters.core_m[slct]), bins=100)
    xbins_cen = dtk.bins_avg(xbins)
    plt.plot(xbins_cen, h, label='>M_infall')


    slct = (clusters.core_m>fit_mi) & (clusters.core_r<fit_rd)
    h, xbins = np.histogram(np.log10(clusters.core_host_mass[slct]/clusters.core_m[slct]), bins=100)
    xbins_cen = dtk.bins_avg(xbins)
    plt.plot(xbins_cen, h, label='Fit Cores')
    
    plt.legend(loc='best')

    
def get_hubble_time(z):
    a = 1.0/(1.0+z)
    H0 = 70*(0.30/a**3 + 0.7) #km/s /MPC
    km_mpc = 3.086e+19 #km/mpc
    h_t = km_mpc/H0
    return h_t/3.154e+7/1e9

def test_wetzel09():
    X = np.logspace(0, 3, 100)
    t_dyn = X/np.log(1.0+X)*0.25
    plt.figure()
    plt.loglog(X, t_dyn, '-')
    

if __name__ == "__main__":

    print("{:.2f}".format(get_hubble_time(0)))
    print("{:.2f}".format(get_hubble_time(.2)))
    print("{:.2f}".format(get_hubble_time(1)))
    # test_wetzel09()
    # calc_wetzel09_disruption(sys.argv[1])


    # calc_disruption(sys.argv[1])
    calc_ngal(sys.argv[1])
    # calc_zeng07()
    plt.show()
