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
import numpy as np



def get_mass_bins(fname):
    hfile = h5py.File(fname,'r')
    return hfile['m_bins'].value

def load_clusters(fname, mass_bins):
    hfile = h5py.File(fname, 'r')
    mass = hfile['cluster/sod_mass'].value
    counts, _ = np.histogram(mass, bins=mass_bins)
    return counts, None


def load_zmrs(fname, mass_bins):
    hfile = h5py.File(fname, 'r')
    mass_bins2 = hfile['m_bins'].value
    assert  np.sum(~(mass_bins == mass_bins2)) == 0, "mass bins are off"
    print(hfile['zm_counts'].value.astype('int'))
    return
    return mass_bin_counts, avg_ngal

def plot_halo_counts():
    OR_clusters_m = "tmp_hdf5/clusters_OR_M200m.hdf5"
    OR_clusters_c = "tmp_hdf5/clusters_OR_M200c.hdf5"
    QC_clusters_m = "tmp_hdf5/clusters_QC_M200m.hdf5"
    QC_clusters_c = "tmp_hdf5/clusters_QC_M200c.hdf5"

    RM_m = "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_simet_mean3/result/type1_weight1_mag1_clr1_result.hdf5"
    RM_c = "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_simet_crit4/result/type1_weight1_mag1_clr1_result.hdf5"
    SP_m = "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_spider/result/type1_weight1_mag1_clr1_result.hdf5"
    SP_c = "/home/dkorytov/phys/Ngal_sdss/data/rad_profile_mstar0_wmap7_spider/result/type1_weight1_mag1_clr1_result.hdf5"

    mass_bins = get_mass_bins(RM_m)

    print("OR Clusters Mean")
    print(load_clusters(OR_clusters_m, mass_bins)[0])
    print("OR Clusters Crit")
    print(load_clusters(OR_clusters_c, mass_bins)[0])
    print("QC Clusters Mean")
    print(load_clusters(QC_clusters_m, mass_bins)[0])
    print("QC Clusters Crit")
    print(load_clusters(QC_clusters_c, mass_bins)[0])

    
    print("Redmapper mean")
    print(load_zmrs(RM_m, mass_bins))
    print("Redmapper crit")
    print(load_zmrs(RM_c, mass_bins))
    print("SPIDERS mean")
    print(load_zmrs(SP_m, mass_bins))
    print("SPIDERS crit")
    print(load_zmrs(SP_c, mass_bins))
          
          




if __name__ == "__main__":
    plot_halo_counts()
