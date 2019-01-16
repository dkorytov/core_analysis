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

        self.core_offset = hfile['cluster/core_offset'].value
        self.core_size = hfile['cluster/core_size'].value
        self.core_x = hfile['cores/core_x'].value
        self.core_y = hfile['cores/core_y'].value
        self.core_z = hfile['cores/core_z'].value
        self.core_r = hfile['cores/core_r'].value
        self.core_m = hfile['cores/core_m'].value

    def plot_cluster(self, i):
        f, (ax1, ax2) = plt.subplots(1,2,figsize=(25,10))
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        ax1.plot(self.core_x[start:stop], self.core_y[start:stop], 'go',mfc='none',mec='g')
        circle = Circle((self.x[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax1.add_artist(circle)
        ax1.grid()
        ax1.set_title("cluster[{}]".format(i))
        ax2.plot(self.core_z[start:stop], self.core_y[start:stop], 'go',mfc='none',mec='g')
        circle = Circle((self.z[i],self.y[i]), self.rad[i], fc='none', ec='k')
        ax2.add_artist(circle)
        ax2.grid()
        ax2.set_title("cores[{}]".format(self.core_size[i]))
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')

    def plot_fit_cluster(self, i , m_infall,  r_disrupt, r_merger = None):
        if self.mass[i] > 2e14:
            return
        start = self.core_offset[i]
        stop = start + self.core_size[i]
        core_x, core_y, core_z  = self.core_x[start:stop], self.core_y[start:stop], self.core_z[start:stop]
        core_r = self.core_r[start:stop]
        core_m = self.core_m[start:stop]
        slct = core_m > m_infall
        slct_compact = core_r < r_disrupt
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
        gal_x, gal_y, gal_z, gal_w = n2merger.n2merger3d(x, y, z, 0.1)
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
        ax1.set_title("cluster[{}] \n Mass {:.2e}".format(i, self.mass[i]))

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
        cluster_data.plot_fit_cluster(i, 12.5, 0.05)
        plt.show()

if __name__ == "__main__":
    plot_saved_clusters(sys.argv[1])



