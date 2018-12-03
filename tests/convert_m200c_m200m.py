#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.cosmology import WMAP7
import halotools.empirical_models




z = np.linspace(0,1, 16)
b = WMAP7.critical_density(z)*10**42
plt.figure()
plt.plot(z, b)
plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200m'), label='200m')
plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200c'), label='200c')
plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '83c'), label='83c')
plt.yscale('log')
plt.legend()


plt.figure()
plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200m')/halotools.empirical_models.density_threshold(WMAP7, z, '200c'))

m200m_to_m200c_density_ratio = halotools.empirical_models.density_threshold(WMAP7, 0.24, '200m')/halotools.empirical_models.density_threshold(WMAP7, 0.24, '200c')

print(m200m_to_m200c_density_ratio)
print(m200m_to_m200c_density_ratio*200)



plt.show()
