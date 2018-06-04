#!/usr/bin/env python2.7

import numpy as np
from catalog_reader import Catalog,frag_to_real
from core_fit2_util import *
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize
from scipy import stats

#plt.rc('font', family='san-serif')
#plt.rc('text', usetex=True)


print "2nd fitting code"
param = dtk.Param(sys.argv[1])

n2lib_loc = param.get_string("n2lib_loc")


fof_loc = param.get_string("fof_loc")
sod_loc = param.get_string("sod_loc")
core_loc = param.get_string("core_loc")
zmr_loc = param.get_string("zmr_loc")

fof_cat = Catalog(fof_loc)
sod_cat = Catalog(sod_loc)
core_cat= Catalog(core_loc)
zmr_sdss = np.load(zmr_loc)

step1 = param.get_int("step")
steps = param.get_int_list("steps")

z_in = param.get_float("z_in")
z_out = param.get_float("z_out")
num_steps = param.get_int("num_steps")
 
stepz = dtk.StepZ(z_in,z_out,num_steps)

n2merger = N2Merger(n2lib_loc)

fof_cat.add_steps(steps)
sod_cat.add_steps(steps)
core_cat.add_steps(steps)


fof_cat.add_data_name("fof_halo_tag")
fof_cat.add_data_name("fof_halo_mass")
fof_cat.add_data_name("fof_halo_center_x")
fof_cat.add_data_name("fof_halo_center_y")
fof_cat.add_data_name("fof_halo_center_z")

sod_cat.add_data_name("fof_halo_tag")
sod_cat.add_data_name("sod_halo_mass")
sod_cat.add_data_name("sod_halo_radius")

core_cat.add_data_name("fof_halo_tag")
core_cat.add_data_name("x")
core_cat.add_data_name("y")
core_cat.add_data_name("z")
core_cat.add_data_name("radius")
core_cat.add_data_name("infall_mass")
print "reading in files"
fof_cat.read_gio()
sod_cat.read_gio()
core_cat.read_gio()

#from the fit. It really should be from a file
#but I'll do that later.
disrupt_len = 0.355
merg_len    = 0.052
infall_mass = 11.61

print "merging catalogs"
halo_cat = Catalog()
halo_cat.join(fof_cat,sod_cat,join_on='fof_halo_tag')
halo_cat.sort('fof_halo_tag')
core_cat.apply_function('fof_halo_tag',frag_to_real)
clstrs = []
zmr_valid = ZMRIndexValidator(zmr_sdss)


htags = halo_cat[step1]['fof_halo_tag']
#how many cores each halo has
fof_cnt = np.zeros(htags.size)
sod_cnt = np.zeros(htags.size)
#selecting the cores have are above the mass cut and that
#are compact enough
core_intct_slct = core_cat[step1]['radius']<disrupt_len
core_infall_slct = core_cat[step1]['infall_mass']>10.0**11.61
core_slct = core_intct_slct & core_infall_slct 


mass_bins = np.logspace(10,16,100)
mass_bins_avg = (mass_bins[:-1]+mass_bins[1:])/2.0
for i in range(0,htags.size):
    if(i%1000 == 0):
        print i,'/',htags.size
    htag = htags[i]
    slct1 = core_cat[step1]['fof_halo_tag']==htag
    slct = slct1 & core_slct
    #TODO account for the mergers
    core_num = np.sum(slct)
    core_x, core_y, core_z,core_w= n2merger.n2merger3d(core_cat[step1]['x'][slct],
                                                   core_cat[step1]['y'][slct],
                                                   core_cat[step1]['z'][slct],
                                                   merg_len)
    fof_cnt[i]=core_xx.size
    #only counting cores within the r200 radius
    sod_x = halo_cat[step1]['fof_halo_center_x'][i]
    sod_y = halo_cat[step1]['fof_halo_center_y'][i]
    sod_z = halo_cat[step1]['fof_halo_center_z'][i]
    sod_r = halo_cat[step1]['sod_halo_radius'][i]
    core_x = core_cat[step1]['x'][slct]
    core_x = core_cat[step1]['y'][slct]
    core_x = core_cat[step1]['z'][slct]
    dx = core_x-sod_x
    dy = core_y-sod_y
    dz = core_z-sod_z
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    within_r200 = dr<sod_r
    core_x, core_y, core_z,core_w
    sod_cnt[i]=within_r200

fof_mean,_,_ = stats.binned_statistic(halo_cat[step1]['fof_halo_mass'],fof_cnt,statistic='mean',bins=mass_bins)
fof_var,_,_  = stats.binned_statistic(halo_cat[step1]['fof_halo_mass'],fof_cnt,statistic=np.std,bins=mass_bins)
fof_cnt,_,_  = stats.binned_statistic(halo_cat[step1]['fof_halo_mass'],fof_cnt,statistic='count',bins=mass_bins)
plt.figure()
plt.plot(mass_bins_avg,fof_mean)
plt.fill_between(mass_bins_avg,fof_mean-fof_var,fof_mean+fof_var,alpha=0.5)
plt.xscale('log')
plt.yscale('log')
plt.title('FoF HOD')
plt.xlabel('FoF halo mass [Msun/h]')
plt.ylabel('Core count')
plt.grid()



sod_mean,_,_ = stats.binned_statistic(halo_cat[step1]['fof_halo_mass'],sod_cnt,statistic='mean',bins=mass_bins)
sod_var,_,_  = stats.binned_statistic(halo_cat[step1]['fof_halo_mass'],sod_cnt,statistic=np.std,bins=mass_bins)
sod_cnt,_,_  = stats.binned_statistic(halo_cat[step1]['sod_halo_mass'],sod_cnt,statistic='count',bins=mass_bins)
plt.figure()
plt.plot(mass_bins_avg,sod_mean)
plt.fill_between(mass_bins_avg,sod_mean-sod_var,sod_mean+sod_var,alpha=0.5)
plt.xscale('log')
plt.yscale('log')
plt.title('SOD HOD')
plt.xlabel('SOD halo mass [Msun/h]')
plt.ylabel('Core count')
plt.grid()

dtk.save_figs("figs/"+param.file+"/plot_hod.py/")

plt.show()
