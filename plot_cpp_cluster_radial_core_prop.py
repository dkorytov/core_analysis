#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy.random as rnd
from scipy import stats
import sys
import h5py
import dtk 


output_loc = "tmp_hdf5/clusters.hdf5"

hfile = h5py.File(output_loc,'r')
cluster_num = hfile["cluster_num"][0]
core_m_host_list = []
core_radial_list = []
core_m_list = []
core_r_list = []
projections = 40
print cluster_num
for i in range(0,cluster_num):
    base_name = "cluster"+str(i)+"/"
    x =  hfile[base_name+"x"][0]
    y =  hfile[base_name+"y"][0]
    z =  hfile[base_name+"z"][0]
    r200 =   hfile[base_name+"sod_radius"][0]
    m200 =   hfile[base_name+"sod_mass"][0]
    core_x = hfile[base_name+"core_x"][:]
    core_y = hfile[base_name+"core_y"][:]
    core_z = hfile[base_name+"core_z"][:]
    core_dr =np.sqrt((core_x-x)**2 + (core_y-y)**2 + (core_z-z)**2)
    #core_dr =np.sqrt((core_x-x)**2 + (core_y-y)**2)
    rand = rnd.rand(projections)*np.pi/2.0
    for i in range(0,projections):
        core_dr200=core_dr/r200*np.cos(rand[i])
        core_m = hfile[base_name+"core_m"][:]
        core_r = hfile[base_name+"core_r"][:]
        core_m_host_list.append(m200*np.ones_like(core_x))
        core_radial_list.append(core_dr200)
        core_m_list.append(core_m)
        core_r_list.append(core_r)


core_m_host = np.concatenate(core_m_host_list)
core_radial = np.concatenate(core_radial_list)
core_m =      np.concatenate(core_m_list)
core_r =      np.concatenate(core_r_list)
print core_radial
m_bins = np.logspace(10,14,100)
m_bins_avg = (m_bins[:-1]+m_bins[1:])/2.0
r_bins = np.logspace(-3,0,100)
r_bins_avg = (r_bins[:-1]+r_bins[1:])/2.0
radial_bins = np.linspace(0,2,25)
radial_bins_avg = (radial_bins[:-1]+radial_bins[1:])/2.0
radial_bins_area = (radial_bins[1:]**2 - radial_bins[:-1]**2)
radial_bins_area = radial_bins_area/radial_bins[-1]**2
m_bins_small = np.logspace(10,14,50)
r_bins_small = np.logspace(-3,0,50)

# plt.figure()
# H,x_bins = np.histogram(core_m_host,bins=np.logspace(11,15,100))
# x_bins_avg = (x_bins[:-1]+x_bins[1:])/2.0
# plt.plot(x_bins_avg,H,'-x')
# plt.yscale('log')
# plt.xscale('log')
# plt.show()


plt.figure()
H,_,_ = np.histogram2d(core_m,core_r,bins=(m_bins,r_bins))
plt.pcolormesh(m_bins,r_bins,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
plt.xscale('log')
plt.yscale('log')
plt.ylabel('core radius [Mpc/h]')
plt.xlabel('central infall mass [Msun/h]')
plt.grid()



plt.figure()
H,_,_ = np.histogram2d(core_radial,core_m,bins=(radial_bins,m_bins))
plt.pcolormesh(radial_bins,m_bins,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial,radial_bins,core_m,np.percentile,16)
median = dtk.binned_func(core_radial,radial_bins,core_m,np.median)
perct84 = dtk.binned_func(core_radial,radial_bins,core_m,np.percentile,84)
plt.plot(radial_bins_avg,perct16,'--r',lw=2)
plt.plot(radial_bins_avg,median,'r',lw=2,label='median')
plt.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
plt.yscale('log')
plt.xlabel('r/r200')
plt.ylabel('Central Infall Mass')
plt.legend(loc='best')
plt.grid()

plt.figure()
H,_,_ = np.histogram2d(core_radial,core_r,bins=(radial_bins,r_bins))
plt.pcolormesh(radial_bins,r_bins,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial,radial_bins,core_r,np.percentile,16)
median = dtk.binned_func(core_radial,radial_bins,core_r,np.median)
perct84 = dtk.binned_func(core_radial,radial_bins,core_r,np.percentile,84)
plt.plot(radial_bins_avg,perct16,'--r',lw=2)
plt.plot(radial_bins_avg,median,'r',lw=2,label='median')
plt.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
plt.yscale('log')
plt.xlabel('r/r200')
plt.ylabel('Core Radius')
plt.legend(loc='best')
plt.grid()




##################################################################################
# Plotting average core prop as function of radial distance from halo center in  #
# bins of host halo mass                                                         #
##################################################################################
labels = []
p16  = []
p50  = []
p84  = []
f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True,figsize=(10,8))

ax = ax1
min_mass = 1e14
max_mass = 2.5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_m[slct],bins=(radial_bins,m_bins_small))
ax.pcolormesh(radial_bins,m_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_ylabel('central infall mass [Msun/h]')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax2
min_mass = 2.5e14
max_mass = 5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_m[slct],bins=(radial_bins,m_bins_small))
ax.pcolormesh(radial_bins,m_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax3
min_mass = 5e14
max_mass = 7.5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_m[slct],bins=(radial_bins,m_bins_small))
ax.pcolormesh(radial_bins,m_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_ylabel('central infall mass [Msun/h]')
ax.set_xlabel('r/r200')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax4
min_mass = 7.5e14
max_mass = 10e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_m[slct],bins=(radial_bins,m_bins_small))
ax.pcolormesh(radial_bins,m_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_m[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel('r/r200')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

clrs = ['b','r','g','c','m']
plt.figure()
for i in range(0,len(labels)):
    plt.plot(radial_bins_avg,p50[i],clrs[i],label=labels[i])
    plt.fill_between(radial_bins_avg,p16[i],p84[i],facecolor=clrs[i],alpha=0.1)
    plt.plot(radial_bins_avg,p16[i],clrs[i]+'--')
    plt.plot(radial_bins_avg,p84[i],clrs[i]+'--')
plt.yscale('log')
plt.legend(loc='best')
plt.ylabel('median central infall mass [Msun/h]')
plt.xlabel('r/r200')

##################################################################################
# Plotting average core radius as function of radial distance from halo center in  #
# bins of host halo mass                                                         #
##################################################################################
labels = []
p16  = []
p50  = []
p84  = []
f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True,figsize=(10,8))

ax = ax1
min_mass = 1e14
max_mass = 2.5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_r[slct],bins=(radial_bins,r_bins_small))
ax.pcolormesh(radial_bins,r_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_ylabel('core radius [Mpc/h]')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax2
min_mass = 2.5e14
max_mass = 5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_r[slct],bins=(radial_bins,r_bins_small))
ax.pcolormesh(radial_bins,r_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()


labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax3
min_mass = 5e14
max_mass = 7.5e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_r[slct],bins=(radial_bins,r_bins_small))
ax.pcolormesh(radial_bins,r_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel('r/r200')
ax.set_ylabel('core radius [Mpc/h]')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

ax = ax4
min_mass = 7.5e14
max_mass = 10e14
slct = (core_m_host >= min_mass) & (core_m_host < max_mass)
ax.set_title('%.2e < Mfof < %.2e'%(min_mass,max_mass))
H,_,_ = np.histogram2d(core_radial[slct],core_r[slct],bins=(radial_bins,r_bins_small))
ax.pcolormesh(radial_bins,r_bins_small,H.T+1,cmap='PuBu',norm=clr.LogNorm())
perct16 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,16)
median = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.median)
perct84 = dtk.binned_func(core_radial[slct],radial_bins,core_r[slct],np.percentile,84)
ax.plot(radial_bins_avg,perct16,'--r',lw=2)
ax.plot(radial_bins_avg,median,'r',lw=2,label='median')
ax.plot(radial_bins_avg,perct84,'--r',lw=2,label='68%')
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
ax.set_xlabel('r/r200')
labels.append('%.2e < Mfof < %.2e'%(min_mass,max_mass))
p16.append(perct16)
p50.append(median)
p84.append(perct84)

clrs = ['b','r','g','c','m']
plt.figure()
for i in range(0,len(labels)):
    plt.plot(radial_bins_avg,p50[i],clrs[i],label=labels[i])
    plt.fill_between(radial_bins_avg,p16[i],p84[i],facecolor=clrs[i],alpha=0.1)
    plt.plot(radial_bins_avg,p16[i],clrs[i]+'--')
    plt.plot(radial_bins_avg,p84[i],clrs[i]+'--')

#plt.yscale('log')
plt.legend(loc='best')
plt.ylabel('median core radius [Mpc/h]')
plt.xlabel('r/r200')



###########################################################################
# Plotting 2d histogram as a function of radial distance from halo center #
###########################################################################
f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True,figsize=(10,8))

rad_min = 0
rad_max = 0.25
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax1.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax1.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid()
ax1.set_ylabel('core radius [Mpc/h]')

rad_min = 0.25
rad_max = 0.50
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax2.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax2.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.grid()

rad_min = 0.50
rad_max = 0.75
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax3.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax3.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.grid()
ax3.set_ylabel('core radius [Mpc/h]')
ax3.set_xlabel('central infall mass [Msun/h]')

rad_min = 0.75
rad_max = 1.0
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax4.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax4.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.grid()
ax4.set_xlabel('central infall mass [Msun/h]')




f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True,figsize=(10,8))

rad_min = 1.00
rad_max = 1.25
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax1.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax1.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid()
ax1.set_ylabel('core radius [Mpc/h]')

rad_min = 1.25
rad_max = 1.50
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax2.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax2.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.grid()

rad_min = 1.50
rad_max = 1.75
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax3.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax3.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.grid()
ax3.set_ylabel('core radius [Mpc/h]')
ax3.set_xlabel('central infall mass [Msun/h]')

rad_min = 1.75
rad_max = 2.0
slct = (core_radial >= rad_min) & (core_radial<rad_max)
ax4.set_title("%.2f< r/r200 < %.2f"%(rad_min,rad_max))
H,_,_ = np.histogram2d(core_m[slct],core_r[slct],bins=(m_bins_small,r_bins_small))
ax4.pcolormesh(m_bins_small,r_bins_small,H.T+1.0,cmap='PuBu',norm=clr.LogNorm())
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.grid()
ax4.set_xlabel('central infall mass [Msun/h]')



###########################
# The abundance of cores  #
###########################

mass_mins = [1e14,2.5e14,5e14,7.5e14]
mass_maxs = [2.5e14,5e14,7.5e14,10e14]

plt.figure()
for i in range(0,len(mass_mins)):
    label = "%.1e < Mfof < %.1e"%(mass_mins[i],mass_maxs[i])
    slct = (core_m_host >= mass_mins[i]) & (core_m_host < mass_maxs[i])
    H, _ = np.histogram(core_radial[slct],bins=radial_bins,normed=True)
    plt.plot(radial_bins_avg,H,label=label)
plt.legend(loc='best')
plt.xlabel('r/r200')
plt.ylabel('core normalized density')



###############################
# The effect on profile shape #
###############################

plt.figure()
m_infall_cuts = [1e11,2e11,5e11,1e12,1e13]
min_mass = 1e14
max_mass = 2.5e14
slct_m = (core_m_host >= min_mass ) & (core_m_host < max_mass)
for i in range(0,len(m_infall_cuts)):
    slct_minfall = core_m > m_infall_cuts[i]
    slct = slct_m & slct_minfall
    num = np.sum(slct)
    H,_ = np.histogram(core_radial[slct],bins = radial_bins,normed=True)
    plt.plot(radial_bins_avg,H/radial_bins_area,label='Mcut=%.1e fract=[%.3f]'%(m_infall_cuts[i],float(num)/float(np.sum(slct_m))))
plt.title("%.1e < Mfof < %.1e"%(min_mass,max_mass))
plt.legend(loc='best')
plt.yscale('log')
plt.xlabel('r/r200')
plt.ylabel('normalized core surface density')

plt.figure()
m_infall_cut = 2e11
r_disrupt_cuts = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1,0.2]
for i in range(0,len(r_disrupt_cuts)):
    slct_rdist = core_r < r_disrupt_cuts[i]
    slct_minfall = core_m < m_infall_cut
    slct = slct_m & slct_rdist & slct_minfall
    num = float(np.sum(slct))
    tot = float(np.sum(slct_m & slct_minfall))
    H,_ = np.histogram(core_radial[slct],bins=radial_bins,normed=True)
    if(np.sum(H) > 0):
        plt.plot(radial_bins_avg,H/radial_bins_area,label='R_disrupt=%.3f fract[%.2f]'%(r_disrupt_cuts[i],num/tot))
plt.title("%.1e < Mfof < %.1e\nM_infall=%.1e"%(min_mass,max_mass,m_infall_cut))
plt.legend(loc='best')
plt.yscale('log')
plt.xlabel('r/r200')
plt.ylabel('normalized core surface density')

plt.figure()
r_disrupt_cut = 0.01
slct_r = core_r < r_disrupt_cut
tot = float(np.sum(slct_r))
for i in range(0,len(m_infall_cuts)):
    slct_m = core_m < m_infall_cuts[i]
    slct = slct_r & slct_m
    num = float(np.sum(slct))
    H,_ = np.histogram(core_radial[slct],bins=radial_bins,normed=True)
    if(np.sum(H)>0):
        plt.plot(radial_bins_avg,H/radial_bins_area,label='M_infall=%.e fract[%.2f]'%(m_infall_cuts[i],num/tot))
plt.title("%.1e < Mfof < %.1e\nR_disrupt=%.3f"%(min_mass,max_mass,r_disrupt_cut))
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel('r/r200')
plt.ylabel('normalized core surface density')

dtk.save_figs(path='figs/radial_core_prop/',extension='.png')
plt.show()
