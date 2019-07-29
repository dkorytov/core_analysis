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
from calc_ngal import *






if __name__ == "__main__":
    mass_bins = np.logspace(12,16,81)
    calc_disruption(sys.argv[1], mass_bins)
