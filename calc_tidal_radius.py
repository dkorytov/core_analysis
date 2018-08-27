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
import astropy 
from astropy.cosmology import WMAP7 as cosmo
import halotools
from halotools.empirical_models import NFWProfile
cosmo
_critical_density_func=None
def rho0(z):
    global _rho0_func
    if _critical_density_func is None:
        zs = np.linspace(0,10,1000)
        rho_0 = 2.77536627e11 
        density = cosmo.critical_density(zs)/cosmo.critical_density0*rho_0
        _critical_density_func = np.interp(z,density)




nfwprofile = NFWProfile()

def F_grav(R,M200,c):
    return dtk.NFW_enclosed_mass(R,c)*M200/R**2

def F_tidal(r, R, M200, c):
    return -F_grav(R,M200,c) + F_grav(R-r, M200, c)

plt.figure()
r = np.linspace(0,0.2,1000)
R = 1
M1 = 1e14
M2 = 1e11
plt.plot(r, F_grav(r,M2,10))
plt.plot(r, F_tidal(r, R, M1,3))
plt.yscale('log')
plt.show()
