#!/usr/bin/env python2.7

from __future__ import print_function, division

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import h5py
import halotools.empirical_models
from astropy.cosmology import WMAP7

def load_halo_cat(sod_location, h_scaling):
    print("loading the catalog")
    cat = {}
    dtk.gio_inspect(sod_location)
    cat['fof_halo_tag'] = dtk.gio_read(sod_location, 'fof_halo_tag'); print("1")
    cat['fof_halo_center_x'] = dtk.gio_read(sod_location, 'fof_halo_center_x'); print("1")
    cat['fof_halo_center_y'] = dtk.gio_read(sod_location, 'fof_halo_center_y'); print("1")
    cat['fof_halo_center_z'] = dtk.gio_read(sod_location, 'fof_halo_center_z'); print("1")
    cat['sod_halo_mass']     = dtk.gio_read(sod_location, 'sod_halo_mass')/h_scaling; print("1")
    cat['sod_halo_radius']     = dtk.gio_read(sod_location, 'sod_halo_radius')/h_scaling; print("1")
    cat['sod_halo_cdelta']           = dtk.gio_read(sod_location, 'sod_halo_cdelta'); print("1")
    slct = cat['sod_halo_cdelta'] < 1
    cat['sod_halo_cdelta'][slct] = 5.5
    slct = cat['sod_halo_cdelta'] > 20
    cat['sod_halo_cdelta'][slct] = 20
    print("done loading")
    return cat


def write_halo_cat(sod_output, cat):
    print("trying to write to ", sod_output)
    hfile = h5py.File(sod_output, 'w')
    hfile['fof_halo_tag'] = cat['fof_halo_tag']
    hfile['fof_halo_center_x'] = cat['fof_halo_center_x']
    hfile['fof_halo_center_y'] = cat['fof_halo_center_y']
    hfile['fof_halo_center_z'] = cat['fof_halo_center_z']
    # M_200m, R_200m
    hfile['sod_halo_mass_m200m']   = cat['sod_halo_mass_m200m']
    hfile['sod_halo_radius_r200m'] = cat['sod_halo_radius_r200m']
    hfile['sod_halo_cdelta_200m']  = cat['sod_halo_cdelta_200m']
    # M_200c, R_200c
    hfile['sod_halo_mass_m200c']   = cat['sod_halo_mass_m200c']
    hfile['sod_halo_radius_r200c'] = cat['sod_halo_radius_r200c']
    hfile['sod_halo_cdelta_200c']  = cat['sod_halo_cdelta']
    

    hfile.close()


def get_rho_mean_over_rho_crit(redshift):
    """Returns the mean density in units of critical density

    """
    return halotools.empirical_models.density_threshold(WMAP7, redshift, '200m')/halotools.empirical_models.density_threshold(WMAP7, redshift, '200c')


def convert_halo_m200c_to_m200m( sod_location, sod_output, time_step, redshift , h_scaling):
    rho_mean = get_rho_mean_over_rho_crit(redshift)
    cat = load_halo_cat(sod_location.replace("${step}", str(time_step)), h_scaling)
    nfw_converter = dtk.NFWConverter(lower_limit = 0.0001, upper_limit = 50)
    starting_delta = 200.0 #
    target_delta = starting_delta*rho_mean
    Rs = cat['sod_halo_radius']/cat['sod_halo_cdelta']
    R_200m = nfw_converter.get_target_overdensity_radius(starting_delta, 
                                                         cat['sod_halo_radius'], 
                                                         cat['sod_halo_cdelta'],
                                                         target_delta,)
    cat['sod_halo_radius_r200m'] = R_200m
    cat['sod_halo_mass_m200m'] = 4.0/3.0 * np.pi * (R_200m**3) * halotools.empirical_models.density_threshold(WMAP7, redshift, '200m')
    cat['sod_halo_cdelta_200m'] = R_200m/Rs
    R_200c = cat['sod_halo_radius']
    cat['sod_halo_mass_m200c'] = 4.0/3.0 * np.pi * (R_200c**3) * halotools.empirical_models.density_threshold(WMAP7, redshift, '200c')
    cat['sod_halo_radius_r200c'] = R_200c
    M_200c = cat['sod_halo_mass_m200c']
    M_200m = cat['sod_halo_mass_m200m']

    h, xbins, ybins = np.histogram2d(np.log10(M_200c), np.log10(M_200m), bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.plot([np.min(xbins), np.max(xbins)], [np.min(xbins), np.max(xbins)], '--k')
    plt.ylabel('M_200m')
    plt.xlabel('M_200c')

    h, xbins, ybins = np.histogram2d(np.log10(M_200c), M_200m/M_200c, bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.ylabel('M_200m/M_200c')
    plt.xlabel('M_200c')

    h, xbins = np.histogram(cat['sod_halo_cdelta'], bins = 100)
    plt.figure()
    plt.plot(dtk.bins_avg(xbins), h, '-x')

    plt.show()
    write_halo_cat(sod_output.replace("${step}", str(time_step)), cat)


if __name__ == "__main__":
    param = dtk.Param(sys.argv[1])
    steps = param.get_int_list("steps")
    sod_location = param.get_string("sod_location")
    sod_output   = param.get_string("sod_output")
    input_h = param.get_float("input_h")
    output_h = param.get_float("output_h")
    h_scaling = input_h/output_h
    stepz = dtk.StepZ(sim_name = 'AlphaQ')
    for step in steps:
        redshift = stepz.get_z(step)
        convert_halo_m200c_to_m200m(sod_location, sod_output, step, redshift, h_scaling)
