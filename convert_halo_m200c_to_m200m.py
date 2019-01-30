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
from colossus.halo import mass_adv
from colossus.cosmology import cosmology


def load_halo_cat(sod_location, h_scaling):
    print("loading the catalog")
    print("H_scaling: ", h_scaling)
    print("mass factor: ", 1/h_scaling)
    cat = {}
    dtk.gio_inspect(sod_location)
    cat['fof_halo_tag'] = dtk.gio_read(sod_location, 'fof_halo_tag'); 
    cat['fof_halo_center_x'] = dtk.gio_read(sod_location, 'fof_halo_center_x'); 
    cat['fof_halo_center_y'] = dtk.gio_read(sod_location, 'fof_halo_center_y'); 
    cat['fof_halo_center_z'] = dtk.gio_read(sod_location, 'fof_halo_center_z'); 
    cat['sod_halo_mass']     = dtk.gio_read(sod_location, 'sod_halo_mass');
    cat['sod_halo_radius']   = dtk.gio_read(sod_location, 'sod_halo_radius');
    cat['sod_halo_cdelta']   = dtk.gio_read(sod_location, 'sod_halo_cdelta'); 
    slct = cat['sod_halo_cdelta'] < 1
    cat['sod_halo_cdelta'][slct] = 5.75
    slct = cat['sod_halo_cdelta'] > 20
    cat['sod_halo_cdelta'][slct] = 20
    print("done loading")
    return cat


def write_halo_cat(sod_output, cat, h_scaling):
    """We convert the radii back to h=0.7 to it will work with GIO hacc
    data.
    """
    print("trying to write to ", sod_output)
    hfile = h5py.File(sod_output, 'w')
    hfile['fof_halo_tag'] = cat['fof_halo_tag']
    hfile['fof_halo_center_x'] = cat['fof_halo_center_x']
    hfile['fof_halo_center_y'] = cat['fof_halo_center_y']
    hfile['fof_halo_center_z'] = cat['fof_halo_center_z']
    # M_200m, R_200m
    hfile['sod_halo_mass_m200m']   = cat['sod_halo_mass_m200m']
    hfile['sod_halo_radius_r200m'] = cat['sod_halo_radius_r200m']*h_scaling
    hfile['sod_halo_cdelta_200m']  = cat['sod_halo_cdelta_200m']
    # M_200c, R_200c
    # hfile['sod_halo_mass_m200c']   = cat['sod_halo_mass_m200c']
    # hfile['sod_halo_radius_r200c'] = cat['sod_halo_radius_r200c']*h_scaling
    # hfile['sod_halo_cdelta_200c']  = cat['sod_halo_cdelta']

    hfile.close()


def get_rho_mean_over_rho_crit(redshift):
    """Returns the mean density in units of critical density

    """
    return halotools.empirical_models.density_threshold(WMAP7, redshift, '200m')/halotools.empirical_models.density_threshold(WMAP7, redshift, '200c')


def convert_halo_m200c_to_m200m( sod_location, sod_output, time_step, redshift , h_scaling, write_catalog=False):
    a = 1.0/(1.0+redshift)
    rho_mean = get_rho_mean_over_rho_crit(redshift)
    cat = load_halo_cat(sod_location.replace("${step}", str(time_step)), h_scaling)
    # nfw_converter = dtk.NFWConverter(lower_limit = 0.0001, upper_limit = 50)
    # starting_delta = 200.0 #
    # target_delta = starting_delta*rho_mean
    # Rs = cat['sod_halo_radius']/cat['sod_halo_cdelta']
    # R_200m = nfw_converter.get_target_overdensity_radius(starting_delta, 
    #                                                      cat['sod_halo_radius'], 
    #                                                      cat['sod_halo_cdelta'],
    #                                                      target_delta,)
    # a = 1.0/(1.0+redshift)
    # cat['sod_halo_radius_r200m'] = R_200m*0.75
    # cat['sod_halo_mass_m200m'] = 4.0/3.0 * np.pi * (R_200m**3) * halotools.empirical_models.density_threshold(WMAP7, redshift, '200m') * 0.7**2 / a**3 
    # cat['sod_halo_cdelta_200m'] = R_200m/Rs
    # R_200c = cat['sod_halo_radius']
    # cat['sod_halo_mass_m200c'] = 4.0/3.0 * np.pi * (R_200c**3) * halotools.empirical_models.density_threshold(WMAP7, redshift, '200c') * 0.7**2 / a**3
    # cat['sod_halo_radius_r200c'] = R_200c*0.75
    # cat['sod_halo_cdelta_200c'] = R_200c/Rs
    # M_200c = cat['sod_halo_mass_m200c']
    # M_200m = cat['sod_halo_mass_m200m']

    M200m_col, R200m_col, c200m_col = mass_adv.changeMassDefinitionCModel(cat['sod_halo_mass']/0.7, redshift,
                                                                          "200c","200m", c_model='child18')#cat['sod_cdelta']
    M200c_col, R200c_col, c200c_col = mass_adv.changeMassDefinitionCModel(cat['sod_halo_mass']/0.7, redshift,
                                                                          "200c","200c", c_model='child18')#cat['sod_cdelta']

    cat['sod_halo_radius_r200m'] = R200m_col/1000.0/0.7*a
    cat['sod_halo_mass_m200m']   = M200m_col*0.7
    cat['sod_halo_cdelta_200m']  = c200m_col

    cat['sod_halo_radius_r200c'] = R200c_col/1000.0/0.7*a
    cat['sod_halo_mass_m200c']   = M200c_col*0.7
    cat['sod_halo_cdelta_200c']  = c200c_col

    M_200c = M200c_col*0.7
    R_200c = R200c_col/1000.0/0.7

    M_200m = M200m_col*0.7
    R_200m = R200m_col/1000.0/0.7

    if( write_catalog):
        print("writing catalog")
        write_halo_cat(sod_output.replace("${step}", str(time_step)), cat, h_scaling)
        exit()
    else:
        print("plotting tests")

    h, xbins, ybins = np.histogram2d(np.log10(cat['sod_halo_mass']), np.log10(M_200c*0.7), bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.plot([np.min(xbins), np.max(xbins)], [np.min(xbins), np.max(xbins)], '--k')
    plt.ylabel('Col M200c')
    plt.xlabel('SOD M200c')

    h, xbins, ybins = np.histogram2d(cat['sod_halo_radius'],R_200m, bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.plot([np.min(xbins), np.max(xbins)], [np.min(xbins), np.max(xbins)], '--k')
    plt.ylabel('Col R200c')
    plt.xlabel('SOD R200c')


    plt.show()
    exit()
    h, xbins, ybins = np.histogram2d(np.log10(M_200c), R_200m/R_200c, bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.ylabel('R_200m/R_200c')
    plt.xlabel('M_200c')

    h, xbins, ybins = np.histogram2d(np.log10(M_200c), M_200m/M_200c, bins = 250)
    plt.figure()
    plt.pcolor(xbins, ybins, h.T, cmap = 'Blues', norm = clr.LogNorm())
    plt.ylabel('M_200m/M_200c')
    plt.xlabel('M_200c')

    h, xbins = np.histogram(cat['sod_halo_cdelta'], bins = 100)
    plt.figure()
    plt.plot(dtk.bins_avg(xbins), h, '-x')

    
    M200m_col, R200m_col, c200m_col = mass_adv.changeMassDefinitionCModel(cat['sod_halo_mass']/0.7, redshift,
                                                                          "200c","200m", c_model='child18')#cat['sod_cdelta']

    plt.figure()
    h, xbins,ybins = np.histogram2d(np.log10(M_200m), np.log10(M200m_col), bins =256)
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    plt.plot([np.min(xbins), np.max(xbins)], [np.min(xbins), np.max(xbins)], '--k')
    plt.xlabel("my m200m")
    plt.ylabel('colossus m200m')


    plt.figure()
    h, xbins, ybins = np.histogram2d(np.log10(M_200c), M_200c/M_200m, bins=100)
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    xbins_cen = dtk.bins_avg(xbins)
    median = dtk.binned_median(np.log10(M_200c), M_200c/M_200m, xbins)
    plt.plot(xbins_cen, median, '-k', label='median')
    plt.legend(loc='best', framealpha=0.0)
    plt.xlabel("M200c")
    plt.ylabel("M200c/M200m")

    plt.figure()
    h, xbins, ybins = np.histogram2d(np.log10(M_200m), M_200c/M_200m, bins=100)
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    xbins_cen = dtk.bins_avg(xbins)
    median = dtk.binned_median(np.log10(M_200m), M_200c/M_200m, xbins)
    plt.plot(xbins_cen, median, '-k', label='median')
    plt.legend(loc='best', framealpha=0.0)
    plt.xlabel("M200m")
    plt.ylabel("M200c/M200m")

    plt.figure()
    M200m_col_lin, R200m_col_lin, c200m_col_lin =  mass_adv.changeMassDefinitionCModel(np.logspace(12,15,100), redshift,
                                                                                       "200c", "200m", c_model='child18')
    h, xbins, ybins = np.histogram2d(np.log10(M_200m), np.log10(R_200m), bins=100)
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    plt.plot(np.log10(M200m_col_lin), np.log10(R200m_col_lin), '--k', label='Colossus')
    plt.xlabel('M200m')
    plt.ylabel('R200m')

    plt.show()


def plot_1to1():
    xlim = plt.xlim()
    ylim = plt.ylim()
    max_val = max(xlim[1], ylim[1])
    min_val = min(xlim[0], ylim[0])
    print(min_val, max_val)
    plt.plot([min_val, max_val], [min_val, max_val], '--k')

def check_halo_conversion(sod_location, time_step, redshift):  
    a = 1.0/(1.0 + redshift)
    cat = load_halo_cat(sod_location.replace("${step}", str(time_step)), h_scaling)
    M200c_col, R200c_col, c200c_col = mass_adv.changeMassDefinitionCModel(cat['sod_halo_mass']/0.7, redshift,
                                                                          "200c","200c", c_model='child18')#cat['sod_cdelta']

    plt.figure()
    plt.plot(cat['sod_halo_mass']/0.7, M200c_col, '.', alpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.title('mass')
    plot_1to1()

    plt.figure()
    plt.plot(cat['sod_halo_radius'], R200c_col/1000.0/0.7*a, '.', alpha=0.3)
    plt.yscale('log')
    plt.xscale('log')
    plt.title('radius')
    plot_1to1()

    plt.figure()
    plt.plot(cat['sod_halo_radius'], (R200c_col/1000.0/0.7*a)/cat['sod_halo_radius'], '.', alpha=0.3)
    plt.xscale('log')
    plt.show()

if __name__ == "__main__":
    cosmology.setCosmology('WMAP7')
    param = dtk.Param(sys.argv[1])
    steps = param.get_int_list("steps")
    sod_location = param.get_string("sod_location")
    sod_output   = param.get_string("sod_output")
    input_h = param.get_float("input_h")
    output_h = param.get_float("output_h")
    write_catalog = param.get_bool("write_catalog")
    h_scaling = input_h/output_h
    stepz = dtk.StepZ(sim_name = 'AlphaQ')
    for step in steps:
        redshift = stepz.get_z(step)
        convert_halo_m200c_to_m200m(sod_location, sod_output, step,
                                    redshift, h_scaling,
                                    write_catalog=write_catalog)
        # check_halo_conversion(sod_location, step, redshift)
