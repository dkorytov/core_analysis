#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
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

class ClusterData:
    def load_file(self, filename):
        hfile = h5py.File(filename, 'r')
        self.num = hfile['cluster_num'].value
        self.x = hfile['cluster/x'].value
        self.y = hfile['cluster/y'].value
        self.z = hfile['cluster/z'].value
        self.rad = hfile['cluster/sod_radius'].value
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
        # ax1.set_aspect('equal')
        # ax2.set_aspect('equal')

        

def plot_saved_clusters(param_filename):
    param = dtk.Param(param_filename)
    cluster_loc = param.get_string("cluster_loc")
    cluster_data = ClusterData();
    cluster_data.load_file(cluster_loc)
    for i in range(0,cluster_data.num):
        cluster_data.plot_cluster(i)
        plt.show()

if __name__ == "__main__":
    plot_saved_clusters(sys.argv[1])


