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

from scipy.special import gammainc


## All in i-band

def press_schecter(norm, Mstar, alpha):
    return lambda M: 0.4*np.log(10) * norm * 10**(-0.4 * (M - Mstar)*(alpha+1.0)) * np.exp(-10**(-0.4*(M-Mstar)))



def luminosity_function_blanton_funct():
    # https://arxiv.org/pdf/astro-ph/0210215.pdf
    return press_schecter(0.0147, -20.82, -1.0)

# def luminosity_function_montero_dorta():
#     return press_schecter(0.0109, -20.97, -1.16)
luminosity_function_montero_dorta = lambda x: press_schecter(0.0109, -20.97, -1.16)(x)
def luminosity_function_blanton(mag_i):
    return luminosity_function_blanton_funct()(mag_i)

if __name__ == "__main__":    
    mags = np.linspace(-24, -15, 1000)
    dmags = mags[1]-mags[0]
    vals = luminosity_function_blanton(mags)
    plt.figure()
    plt.plot(mags, np.log10(luminosity_function_montero_dorta(mags)))
    plt.grid()

    lum_func = luminosity_function_montero_dorta(mags)
    cum_sum = np.cumsum(lum_func)*dmags
    plt.figure()
    plt.plot(mags, cum_sum)
    plt.yscale('log')
    
    interpolation = np.interp(-21.22, mags, cum_sum)
    print(interpolation)
    plt.show()

    
    
