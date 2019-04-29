#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import os
import sys
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr

import dtk
import h5py
from matplotlib.patches import Circle
from core_fit2_util import *


class ClusterData:

    def load_file(self, filename):
        hfile = h5py.File(filename, 'r')
        self.num = hfile['cluster_num'].value
        self.x =    hfile['cluster/x'].value
        self.y =    hfile['cluster/y'].value
        self.z =    hfile['cluster/z'].value
        self.rad =  hfile['cluster/sod_radius'].value
        self.mass = hfile['cluster/sod_mass'].value
        self.htag = hfile['cluster/htag'].value

        self.core_offset = hfile['cluster/core_offset'].value
        self.core_size = hfile['cluster/core_size'].value
        self.core_htag = hfile['cores/core_htag'].value
        self.core_infall_htag = hfile['cores/core_infall_htag'].value
        print(hfile['cores'].keys())
        self.core_x = hfile['cores/core_x'].value
        self.core_y = hfile['cores/core_y'].value
        self.core_z = hfile['cores/core_z'].value
        self.core_r = hfile['cores/core_r'].value
        self.core_m = hfile['cores/core_m'].value
        self.core_step = hfile['cores/core_step'].value
        #self.core_infall_htag= hfile['cores/infall_htag'].value
        self.core_is_central = hfile['cores/core_is_central'].value
        

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
        plt.close()

    def plot_find_central(self, i):
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        slct = core_m > m_infall
        slct_compact = core_r < r_disrupt
        
    def plot_fit_cluster(self, i , m_infall,  r_disrupt, r_merger = 0.0):
        if self.mass[i] > 2e14:
            return
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        slct = core_m > m_infall
        slct_compact = core_r < r_disrupt
        
        dx = core_x - self.x[i]
        dy = core_y - self.y[i]
        dz = core_y - self.z[i]
        dr = dx*dx + dy*dy + dz*dz
        
        srt = np.argsort(dr)
        print("{:.2e}".format(np.max(core_m[srt][:100])))
        # f, (ax1, ax2) = plt.subplots(1,2,figsize=(25,10))
        # ax1.plot(core_x[slct & ~slct_compact], core_y[slct & ~slct_compact], 'go',mfc='none',mec='g')
        # ax1.plot(core_x[slct &  slct_compact], core_y[slct &  slct_compact], 'go',mfc='g',mec='g')
        # circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        # ax1.add_artist(circle)
        # ax1.grid()
        # 
        # ax2.plot(core_z[slct & ~slct_compact], core_y[slct & ~slct_compact], 'go',mfc='none',mec='g')
        # ax2.plot(core_z[slct &  slct_compact], core_y[slct &  slct_compact], 'go',mfc='g',mec='g')
        # circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        # ax2.add_artist(circle)
        # ax2.grid()
        # 
        # ax1.set_aspect('equal')
        # ax2.set_aspect('equal')
        # f = plt.figure()
        # ax = f.gca(projection = "3d")
        # ax.scatter(core_x[slct &  slct_compact], core_y[slct &  slct_compact], core_z[slct &  slct_compact], edgecolors='r', facecolor = 'r', s = 3, alpha=0.3)
        # ax.scatter(core_x[slct & ~slct_compact], core_y[slct & ~slct_compact], core_z[slct & ~slct_compact], edgecolors='g', facecolor = 'none', s=3, alpha=0.3)
        # u = np.linspace(0, 2 * np.pi, 100)
        # v = np.linspace(0, np.pi, 100)

        # x = self.rad[i] * np.outer(np.cos(u), np.sin(v)) + self.x[i]
        # y = self.rad[i] * np.outer(np.sin(u), np.sin(v)) + self.y[i]
        # z = self.rad[i] * np.outer(np.ones(np.size(u)), np.cos(v)) + self.z[i]
        
        # ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.2)
        # # Create cubic bounding box to simulate equal aspect ratio
        # max_range = np.array([core_x.max()-core_x.min(), core_y.max()-core_y.min(), core_z.max()-core_z.min()]).max()
        # Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(core_x.max()+core_x.min())
        # Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(core_y.max()+core_y.min())
        # Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(core_z.max()+core_z.min())
        # # Comment or uncomment following both lines to test the fake bounding box:
        # for xb, yb, zb in zip(Xb, Yb, Zb):
        #     ax.plot([xb], [yb], [zb], 'w')

        # ax.set_aspect('equal')
        
        f, (ax1, ax2) = plt.subplots(1, 2, figsize = (25, 10), sharey=True)
        x = core_x[slct &  slct_compact]
        y = core_y[slct &  slct_compact]
        z = core_z[slct &  slct_compact]
        gal_x, gal_y, gal_z, gal_w = n2merger.n2merger3d(x, y, z, r_merger)
        ax1.scatter(x, y, facecolor='g', edgecolor='g', alpha=0.3)
        ax1.scatter(gal_x, gal_y, facecolor='r', edgecolor='r', alpha=0.3)

        ax2.scatter(z, y, facecolor='g', edgecolor='g', alpha=0.3)
        ax2.scatter(gal_z, gal_y, facecolor='r', edgecolor='r', alpha=0.3)

        circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax1.add_artist(circle)
        circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax2.add_artist(circle)
        
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax2.set_title("cores[{}]".format(self.core_size[i]))
        ax1.set_title("cluster[{}] \n Mass {:.2e}, R200: {:.2f}".format(i, self.mass[i], self.rad[i]))


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
        

n2merger= None

def plot_saved_clusters(param_filename):
    global n2merger
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string("cluster_loc")
    cluster_data = ClusterData();
    n2lib_loc = "lib/libn2merg.so"
    n2merger = N2Merger(n2lib_loc)

    print(n2lib_loc)
    cluster_data.load_file(cluster_loc)
    for i in range(0,cluster_data.num):
        # cluster_data.plot_fit_cluster(i, 12.2, 0.02)
        cluster_data.plot_cluster(i)
        plt.show()

def plot_hist(data, bins, style='-', label=None, ):
    h, xbins = np.histogram(data, bins=bins)
    bins_cen = dtk.bins_avg(xbins)
    plt.plot(bins_cen, h, style, label=label)

def test_saved_clusters(param_filename):
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string('cluster_loc')
    sod_loc = param.get_string("sod_loc")
    sod_hdf5 = param.get_string("sod_hdf5")
    
    sod_data_raw = SODData();

    sod_data_saved = SODData();
    sod_data_saved.load_sod('/media/luna1/dkorytov/data/OuterRim/sod_200m/sod_m200m.401.test1.hdf5', True)

    saved_cluster = ClusterData()
    saved_cluster.load_file('tmp_hdf5/clusters_OR_M200m.test1.hdf5')
    # sod_data_raw.load_sod('/media/luna1/dkorytov/data/OuterRim/sod/m000.401.sodproperties', False)
    plt.figure()
    bins = np.logspace(12, 16, 32)
    plt.title("Mass")
    # plot_hist(sod_data_raw.mass, bins, label='sod_raw')
    plot_hist(sod_data_saved.mass, bins, label='sod saved')
    plot_hist(saved_cluster.mass, bins, label='saved cluster')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('m200')
    plt.legend(loc='best')

    plt.figure()
    bins = np.linspace(0,5,64)
    plt.title("radius")
    # plot_hist(sod_data_raw.rad, bins, label='sod_raw')
    plot_hist(sod_data_saved.rad, bins, label='sod saved')
    plot_hist(saved_cluster.rad, bins, label='saved cluster')
    plt.yscale('log')
    plt.xlabel('r200')
    plt.legend(loc='best')
    plt.show()

if __name__ == "__main__":
    plot_saved_clusters(sys.argv[1])
    #test_saved_clusters(sys.argv[1])


