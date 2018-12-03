#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.cosmology import WMAP7
import halotools.empirical_models
from scipy.interpolate import interp1d



# z = np.linspace(0,1, 16)
# b = WMAP7.critical_density(z)*10**42
# plt.figure()
# plt.plot(z, b)
# plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200m'), label='200m')
# plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200c'), label='200c')
# plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '83c'), label='83c')
# plt.yscale('log')
# plt.legend()


# plt.figure()
# plt.plot(z, halotools.empirical_models.density_threshold(WMAP7, z, '200m')/halotools.empirical_models.density_threshold(WMAP7, z, '200c'))

m200m_to_m200c_density_ratio = halotools.empirical_models.density_threshold(WMAP7, 0.24, '200m')/halotools.empirical_models.density_threshold(WMAP7, 0.24, '200c')

print(m200m_to_m200c_density_ratio)
print(m200m_to_m200c_density_ratio*200)

# def nfw_enclosed_density(m000c, radius, conc):


    

nfw_model = halotools.empirical_models.NFWProfile()

radius = np.linspace(0,2.0, 100)
enclosed_mass = nfw_model.enclosed_mass(radius, 1e14, 3.0)
Rs = nfw_model.halo_mass_to_halo_radius(1e14)
print(Rs, Rs/3.0)
plt.figure()
plt.plot(radius, enclosed_mass)
plt.yscale('log')

plt.figure()
plt.plot(radius, enclosed_mass)
plt.yscale('log')



def nfw_mass_enclosed_scale(r, Rs):
    return np.log((Rs+r)/Rs) - r/(Rs+r)

def nfw_density_enclosed_scale(r, Rs):
    return nfw_mass_enclosed_scale(r, Rs)/r**3

class NFWConverter:
    def __init__(self,):
        radius = np.linspace(0.01, 10, 5000)
        #Unitless density
        enclosed_density = nfw_density_enclosed_scale(radius, 1.0)
        log10_enclosed_density = np.log10(enclosed_density)
        # print("max/min of enclosed density: ",np.max(log10_enclosed_density), np.min(log10_enclosed_density))
        self.enclosed_density_from_radius = interp1d(radius, log10_enclosed_density)
        self.radius_from_enclosed_density = interp1d(log10_enclosed_density, radius)

    def get_target_overdensity_radius(self, starting_delta, R_delta, conc, target_delta):
        """We take a nfw halo with a known R200c and concetration and find
        the R(delta)c, where delta is a specified over density. The function
        return the radius were the target overdensity is reached"""
        # print("we are converting R_{:.0f} = {:.2f}, c={:.2f} to R_{:.0f}".format(starting_delta, R_delta, conc, target_delta))
        R_delta_Rs_units = R_delta*conc
        # print("R_{:.0f} in units of Rs is {:.2f}".format(starting_delta, R_delta_Rs_units))
        starting_density = self.enclosed_density_from_radius(R_delta_Rs_units)
        # print("The average over log10 density at R_{:.0f} is {:.2f}".format(starting_delta, starting_density))
        target_density = starting_density+np.log10(target_delta/starting_delta)
        # print("If we are going from R_{:.0f} to R_{:.0f}, the new log10 density should be {:.2f}".format(starting_delta, target_delta, target_density))
        R_target_delta_Rs_units = self.radius_from_enclosed_density(target_density)
        return R_target_delta_Rs_units/conc
        
                                                      
# plt.figure()
# plt.plot(radius, nfw_mass_enclosed_scale(radius, 0.319175317807)/enclosed_mass)

# plt.figure()
# plt.plot(radius, nfw_mass_enclosed_scale(radius, 0.3))
# plt.plot(radius, nfw_mass_enclosed_scale(radius/0.3, 1.0), '--r')

# plt.figure()
# plt.plot(radius, nfw_density_enclosed_scale(radius, 0.3))
# plt.yscale('log')

nfw_converter = NFWConverter()




print(nfw_converter.get_target_overdensity_radius(200, 1.0, 3.0, 83))




plt.show()
