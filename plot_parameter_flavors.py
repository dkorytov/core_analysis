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


def insert_value(dict, param_name, val, err_p, err_m, mass_def, fit_type, model):
    param = dict.setdefault(mass_def, {}).setdefault(fit_type, {}).setdefault(model, {}).setdefault(param_name, {})

    param['val'] = val
    param['err_p'] = err_p
    param['err_m'] = err_m
    

def create_parameter_record():
    """
    The record goes as mstar->mass_def->fit_type->model_flavor
    """
    result = {}
    insert_value(result,"M_infall", 12.13, .1, -.1, 'simet_mean', 'profile', 'D')
    insert_value(result,"R_disrpt", 0.024, .1, -.1, 'simet_mean', 'profile', 'D')
    return result


def plot_parameter_flavors():
    rec = create_parameter_record()
    print(rec)

if __name__ == "__main__":
    plot_parameter_flavors()

