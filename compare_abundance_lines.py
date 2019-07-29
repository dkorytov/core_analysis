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

import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm

def plot_abundance_line(ax, param_fname):
    param = dtk.Param(param_fname)
    core_loc = param.get_string('core_loc')
    expected_abundance= param.get_string('expected_abundance')
    fname = "tmp_hdf5/{}/abundance={}.hdf5".format(core_loc, expected_abundance)
    print(fname)

def compare_abudance_lines(param_fnames):
    f, ax = plt.subplots(1,1)
    for param_fname in param_fnames:
        plot_abundance_line(ax, param_fname)


if __name__ == "__main__":
    compare_abundance_lines(sys.argv[1:])
