#!/usr/bin/env python3

from __future__ import print_function, division 
import numpy as np
import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import sys
import dtk
import h5py


def load_prop(fname, var, cat):
    cat[var] = dtk.gio_read(fname, var)
    
def load_cores(fname, step=None, quick=False):
    cat = {}
    # dtk.gio_inspect(fname)
    load_prop(fname, 'core_tag', cat)
    load_prop(fname, 'x', cat)
    load_prop(fname, 'y', cat)
    load_prop(fname, 'z', cat)
    load_prop(fname, 'infall_tree_node_mass', cat)
    load_prop(fname, 'infall_step', cat)
    load_prop(fname, 'central', cat)
    if not quick:
        load_prop(fname, 'vx', cat)
        load_prop(fname, 'vy', cat)
        load_prop(fname, 'vz', cat)
        load_prop(fname, 'radius', cat)
        load_prop(fname, 'infall_fof_halo_mass', cat)

    if step is not None:
        size = cat['x'].size
        cat['step'] = np.ones(size)*step
        
    return cat

def get_axis_bins(data, log=False, num=100):
    upper, lower = np.max(data), np.min(data)
    if log:
        upper, lower = int(np.log10(upper))+1, int(np.log10(lower))-1
        return np.logspace(lower, upper, num)
    else:
        return num
def param_set_true(log):
    if log is not None:
        return log
    else:
        return False
        
def plot_2d_dist(cat, var1, var2, bins=None, log_all=False, log1=False, log2=False, central=None):
    plt.figure()
    x = cat[var1]
    y = cat[var2]
    if central is not None:
        if central:
            slct = cat['central']==1
        else:
            slct = cat['central']==0
        x = x[slct]
        y = y[slct]
        plt.title('Central={}'.format(central))
    if log_all:
        log1 = True
        log2 = True
    
    xbins = get_axis_bins(x, log1)
    ybins = get_axis_bins(y, log2)
    h, xbins, ybins = np.histogram2d(x, y, bins=(xbins, ybins))
    plt.pcolor(xbins, ybins, h.T, cmap='Blues', norm=clr.LogNorm())
    if log1:
        plt.xscale('log')
    if log2:
        plt.yscale('log')
    plt.xlabel(var1)
    plt.ylabel(var2)

def plot_core_properties(core_fname):
    cat = load_cores(core_fname)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'infall_step', log1=True)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'radius', log_all=True)
    plot_2d_dist(cat, 'infall_step', 'radius', log_all=True)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'infall_step', log1=True, central=True)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'radius', log_all=True, central=True)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'infall_step', log1=True, central=False)
    plot_2d_dist(cat, 'infall_tree_node_mass', 'radius', log_all=True, central=False)
    plot_2d_dist(cat, 'x', 'y')
    plot_2d_dist(cat, 'x', 'z')
    plot_2d_dist(cat, 'vx', 'vy')
    plot_2d_dist(cat, 'vx', 'vz')
    plot_2d_dist(cat, 'x', 'vx')
    plot_2d_dist(cat, 'x', 'vy')
    plot_2d_dist(cat, 'x', 'vz')

def get_full_core_catalog(core_fname_pattern, steps):
    cats = {}
    for step in steps:
        print(step)
        cat = load_cores(core_fname_pattern.replace('$step$', str(step)), step)
        cats[step] = cat
    final_cat = {}
    keys = cats[steps[0]].keys()
    for key in keys:
        column = []
        for step in steps:
            column.append(cats[step][key])
        final_cat[key]=np.concatenate(column)
    return final_cat
    
def plot_core_trajectory(core_fname_pattern, steps):
    cat = get_full_core_catalog(core_fname_pattern, steps)

    slct = cat['step']==steps[0]
    print('select size: ', slct.size/ np.sum(slct))

    # sort ids by infall tree node mass
    core_masses = cat['infall_tree_node_mass'][slct]
    srt = np.argsort(core_masses)[::-1]
    srt = np.argsort(srt)
    core_tags = cat['core_tag'][slct][srt]
    core_masses = core_masses[srt]

    for core_tag, core_mass in zip(core_tags, core_masses):
        slct = cat['core_tag']==core_tag
        srt = np.argsort(cat['step'][slct])
        x = cat['x'][slct][srt]
        y = cat['y'][slct][srt]
        z = cat['z'][slct][srt]
        vx = cat['vx'][slct][srt]
        vy = cat['vy'][slct][srt]
        vz = cat['vz'][slct][srt]
        step = cat['step'][slct][srt]
        radius = cat['radius'][slct][srt]
        mass = cat['infall_tree_node_mass'][slct][srt]
        central = cat['central'][slct][srt] == 1
        print(list(zip(x,y)))
        print(cat['step'][slct][srt])
        step_mid = dtk.bins_avg(step)
        central_step = np.diff(central)
        
        fig, axs = plt.subplots(2, 3, figsize=(15, 8))
        fig.suptitle('TreeNodeMass={:.2e}\nCore Tag:{}'.format(core_mass, core_tag))
        axs[0][0].plot(x-x[-1], y-y[-1], '-o', mfc='none', label='Path')
        axs[0][0].plot(x[-1]-x[-1], y[-1]-y[-1], 'o', ms=20, mfc='none', label='Last Pos', color='tab:orange')
        axs[0][0].plot(x[central]-x[-1], y[central]-y[-1], 'o', label='Central', color='tab:blue')

        axs[0][0].legend(framealpha=0.0)
        axs[0][0].axis('equal')
        axs[0][0].set_xlabel('x')
        axs[0][0].set_ylabel('y')
        
        axs[0][1].plot(z-z[-1], y-y[-1], '-o', mfc='none', label='Path')
        axs[0][1].plot(z[-1]-z[-1], y[-1]-y[-1], 'o', ms=20, mfc='none', label='Last Pos')
        axs[0][1].plot(z[central]-z[-1], y[central]-y[-1], 'o', label='Central', color='tab:blue')
        axs[0][1].axis('equal')
        axs[0][1].set_xlabel('z')
        axs[0][1].set_ylabel('y')
        axs[0][1].set_ylim(axs[0][0].get_ylim())

        dx = np.diff(x)
        dy = np.diff(y)
        dz = np.diff(z)
        dr = np.sqrt(dx*dx + dy*dy + dz*dz)

        axs[0][2].plot(step_mid, dr, '-o')
        axs[0][2].fill_between(step_mid, 0, dr, alpha=0.3)
        axs[0][2].set_xlabel('Mid Step')
        axs[0][2].set_ylabel('displace per snapshot [Mpc/h]')

        dv = np.sqrt(vx*vx + vy*vy + vz*vz)[1:]
        # dx = vx[1:]
        # dy = vy[1:]
        # dz = vz[1:]
        theta = np.arccos((vx[1:]*dx + vy[1:]*dy + vz[1:]*dz)/(dv*dr))
        theta = theta*180.0/np.pi
        
        # theta = (dx*dx + dy*dy + dz*dz)/(dr*dr)
        # theta = np.arccos
        axs[1][2].plot(step_mid, theta, '-o')
        axs[1][2].fill_between(step_mid, 0, theta, alpha=0.3)
        axs[1][2].set_xlabel('Mid Step')
        axs[1][2].set_ylabel('Angle between displ. and vel. [deg.]')

        axs[1][0].plot(step, mass, '-o', mfc='none')
        axs[1][0].fill_between(step, np.min(mass), mass, alpha=0.3)
        axs[1][0].plot(step[-1], mass[-1], 'o', ms=20, mfc='none')
        axs[1][0].plot(step[central], mass[central], 'o', color='tab:blue')
        axs[1][0].set_xlabel('Step')
        axs[1][0].set_ylabel('Infall Mass [Msun/h]')
        axs[1][0].set_yscale('log')

        axs[1][1].plot(step, radius*1000.0, '-o', mfc='none')
        axs[1][1].fill_between(step, 0, radius*1000.0, alpha=0.3)
        axs[1][1].plot(step[-1], radius[-1]*1000.0, 'o', ms=20, mfc='none')
        axs[1][1].plot(step[central], radius[central]*1000.0, 'o', color='tab:blue')
        axs[1][1].set_xlabel('Step')
        axs[1][1].set_ylabel('Core Radius [Kpc/h]')
        axs[1][1].set_ylim([0, 200])
        # axs[1][1].set_yscale('log')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        # plt.show()
        dtk.save_figs('figs/'+__file__+'/individual_cores_random_rev4070_v2/')
        plt.close()
        



if __name__ == "__main__":
    steps=[499, 487, 475, 464, 453, 442, 432, 421, 411, 401, 392, 382, 373, 365, 355, 347, 338, 331, 323, 315, 307, 300, 293, 286, 279, 272, 266, 259, 253, 247, 241, 235, 230, 224, 219, 213, 208, 203, 198, 194, 189, 184, 180, 176, 171, 163, 159, 155, 151, 148, 144, 141, 137, 134, 131, 127, 124, 121, 119, 116, 113, 110, 107, 105, 102, 100, 97, 95, 92, 90, 88, 84, 81, 79, 77, 76, 74, 72, 70, 68, 67, 65, 63, 62, 60, 59, 57, 56, 54, 53, 52, 50, 49, 48, 46, 45, 44, 43]
    # steps=[307, 300, 293, 286, 279, 272, 266, 259, 253, 247, 241, 235, 230, 224, 219, 213, 208, 203, 198, 194, 189, 184, 180, 176, 171, 163, 159, 155, 151, 148, 144, 141, 137, 134, 131, 127, 124, 121, 119, 116, 113, 110, 107, 105, 102, 100]
    # steps=[307, 300, 293, 286, 279, 272, 266, 259, 253, 247, 241, 235, 230, 224, 219, 213, 208]
    # steps=[307, 300, 293, 286, 279, 272, 266, 259]
    core_fname_pattern = '/data/a/cpac/dkorytov/data/AlphaQ/core_catalog10_rev4070/$step$.coreproperties'
    # plot_core_properties(core_fname_pattern.replace('$step$', str(487)))
    # dtk.save_figs('figs/'+__file__+'/')
    plt.show()
    plot_core_trajectory(core_fname_pattern, steps)

