#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import h5py
import dtk
import sys
import corner
from plot_mcmc_rg import gelman_ruben_test
param = dtk.Param(sys.argv[1])

param_cnt = 2
fit_r_fof= param.get_bool("fit_r_fof")
fit_r_merger= param.get_bool("fit_r_merger")
if(fit_r_fof):
    param_cnt +=1
if(fit_r_merger):
    param_cnt +=1

def autocorr(x):
    x2 = x-np.average(x)
    x3 = x2/np.sqrt(np.sum(x2*x2))
    result = np.correlate(x3,x3,mode='full')
    return (np.arange(0,x.size-1),result[x.size:])

mcmc_loc = "output/"+param.file_name+"/mcmc.gio"
#mcmc_loc = '/home/dkorytov/tmp/mcmc.gio'
#mcmc_loc = "jupiter_result/mcmc.gio"

print "loading mi"
mcmc_m_i  = dtk.gio_read(mcmc_loc,"mcmc_mass_infall")
print "loading rd"
mcmc_r_d  = dtk.gio_read(mcmc_loc,"mcmc_r_disrupt")
print "loading id"
mcmc_id   = dtk.gio_read(mcmc_loc,"mcmc_walker_id")
print "loading ws"
mcmc_step = dtk.gio_read(mcmc_loc,"mcmc_walker_step")
print "loading val"
mcmc_val  = dtk.gio_read(mcmc_loc,"mcmc_value")
hfile = h5py.File("output/"+param.file_name+"/fit_core_params.hdf5")
fitted_m_infall = hfile['m_infall']
fitted_r_disrupt = hfile['r_disrupt']

def stats(data):
    mean = np.average(data)
    std = np.std(data)
    half = data.size/2
    mean2 = np.average(data[half:])
    std2  = np.std(data[half:])
    return mean,std,mean2,std2

def stats2(data,mcmc_id,mcmc_step):
    print data.shape
    print mcmc_id.shape
    print mcmc_step.shape
    max_step = np.max(mcmc_step)
    slct_half = mcmc_step>max_step/2.0
    mean2 = np.average(data[slct_half])
    std2  = np.std(data[slct_half])
    rg2 = gelman_ruben_test(data[slct_half],mcmc_id[slct_half])
    return rg2,mean2,std2

for i in range(0,10):
    print "%d[%d]: %.2e %.3f"%(mcmc_id[i],mcmc_step[i], mcmc_m_i[i], mcmc_r_d[i])



size = np.max(mcmc_step)
slct_half = mcmc_step>size/2
mcmc_mi_lg = np.log10(mcmc_m_i[slct_half])
mcmc_rd_lg = np.log10(mcmc_r_d[slct_half])
H, x_bins, y_bins = np.histogram2d(mcmc_mi_lg,mcmc_rd_lg,bins = (100,100))
H = H.astype(float)
x_bins_avg = (x_bins[0:-1]+x_bins[1:])/2.0
y_bins_avg = (y_bins[0:-1]+y_bins[1:])/2.0
plt.figure(figsize=(10,8))
ax1 = plt.subplot(221)
ax3 = plt.subplot(223,sharex=ax1)
ax4 = plt.subplot(224,sharey=ax3)

ax1.plot(x_bins_avg,np.sum(H,axis=1)/np.sum(H))
ax1.set_xlabel(r'log10(M$_{infall}$)')
ax3.pcolor(x_bins,y_bins,H.T,cmap='PuBu',norm = clr.LogNorm())
ax3.scatter(np.log10(fitted_m_infall),np.log10(fitted_r_disrupt),facecolor='none',edgecolor='r')
ax3.set_xlabel(r'log10(M$_{infall}$)')
ax3.set_ylabel(r'log10(R$_{disrupt}$)')
ax4.plot(np.sum(H,axis=0)/np.sum(H),y_bins_avg)
ax4.grid(True)
ax3.grid(True)
ax1.grid(True)

plt.figure()
max_id = np.max(mcmc_id)
for i in range(0,max_id+1):
    slct = mcmc_id == i
    plt.plot(mcmc_step[slct],mcmc_m_i[slct],'-',alpha=0.3)
plt.ylabel('M_infall [Msun/h]')
plt.xlabel('step number')
plt.yscale('log')
plt.title("chain m_infall \n1/2 rg:%.4f mean:%.4f,std:%.4f"%stats2(np.log10(mcmc_m_i),mcmc_id,mcmc_step))
#plt.show()

plt.figure()
max_id = np.max(mcmc_id)
acpt_rate = []
for i in range(0,max_id+1):
    slct = mcmc_id == i
    print i 
    tmp = mcmc_m_i[slct][::param_cnt]
    x,corr = autocorr(tmp)
    tmp_diff = tmp[1:]-tmp[:-1]
    tmp_zero = tmp_diff==0
    acpt_rate.append(1-np.average(tmp_zero))
    plt.plot(x+1.0,corr)
plt.title('auto corr M_infall')
plt.xscale('log')

#plt.show()                     
print "Average acceptance rate m_i: ", np.average(acpt_rate)


plt.figure()
max_id = np.max(mcmc_id)
for i in range(0,max_id+1):
    slct = mcmc_id == i
    plt.plot(mcmc_step[slct],mcmc_r_d[slct],'-',alpha=0.3)
plt.title("chain r_disrupt \n1/2 rg:%.4f mean:%.4f,std:%.4f"%stats2(np.log10(mcmc_r_d),mcmc_id,mcmc_step))
plt.ylabel('R_disrupt [Mpc/h]')
plt.xlabel('step number')
plt.yscale('log')


plt.figure()
max_id = np.max(mcmc_id)
acpt_rate = []
for i in range(0,max_id+1):
    slct = mcmc_id == i
    tmp = mcmc_r_d[slct][::param_cnt]
    x,corr = autocorr(tmp)
    tmp_diff = tmp[1:]-tmp[:-1]
    tmp_zero = tmp_diff==0
    acpt_rate.append(1-np.average(tmp_zero))
    plt.plot(x+1.0,corr)
plt.title("auto corr R_dirsupt")
plt.xscale('log')
print "Average acceptance rate r_d: ", np.average(acpt_rate)


mi = np.log10(mcmc_m_i)
rd = np.log10(mcmc_r_d)
fit = np.polyfit(mi,rd,1)
print fit

mi = mi-np.average(mi)
mi = mi/np.sqrt(np.sum(mi*mi))

rd = rd-np.average(rd)
rd = rd/np.sqrt(np.sum(rd*rd))
#rd = rd/np.sqrt(np.sum(rd*rd))

result =  np.correlate(mi,mi)

if(param.get_bool("fit_r_fof")):
    mcmc_r_fof = dtk.gio_read(mcmc_loc,"mcmc_r_fof")
    plt.figure()
    max_id = np.max(mcmc_id)
    for i in range(0,max_id+1):
        slct = mcmc_id == i
        plt.plot(mcmc_step[slct],mcmc_r_fof[slct],'-',alpha=0.3)
    plt.ylabel('R_disrupt [Mpc/h]')
    plt.xlabel('step number')

    plt.figure()
    max_id = np.max(mcmc_id)
    acpt_rate = []
    for i in range(0,max_id+1):
        slct = mcmc_id == i
        tmp = mcmc_r_fof[slct]
        x,corr = autocorr(tmp)
        tmp_diff = tmp[1:]-tmp[:-1]
        tmp_zero = tmp_diff==0
        acpt_rate.append(1-np.average(tmp_zero))
        plt.plot(x+1.0,corr)
    plt.title("AutoCorr R_fof")
    print "Average acceptance rate r_fof: ", np.average(acpt_rate)

if(param.get_bool("fit_r_merger")):
    mcmc_r_merger = dtk.gio_read(mcmc_loc,"mcmc_r_merger")
    print "looking at mcmc_r_merger..."
    plt.figure()
    max_id = np.max(mcmc_id)
    for i in range(0,max_id+1):
        slct = mcmc_id == i
        plt.plot(mcmc_step[slct],mcmc_r_merger[slct],'-',alpha=0.3)
    plt.ylabel('R_merger [Mpc/h]')
    plt.yscale('log')
    plt.xlabel('step number')

    plt.figure()
    max_id = np.max(mcmc_id)
    acpt_rate = []
    for i in range(0,max_id+1):
        slct = mcmc_id == i
        tmp = mcmc_r_merger[slct][::param_cnt]
        x,corr = autocorr(tmp)
        tmp_diff = tmp[1:]-tmp[:-1]
        tmp_zero = tmp_diff==0
        acpt_rate.append(1-np.average(tmp_zero))
        plt.plot(x+1.0,corr)
    plt.xscale('log')
    plt.title("AutoCorr R_merger")
    print "Average acceptance rate r_merger: ", np.average(acpt_rate)

    data = np.column_stack((np.log10(mcmc_m_i[slct_half]),np.log10(mcmc_r_d[slct_half]),np.log10(mcmc_r_merger[slct_half])))
    
    figure = corner.corner(data,labels=[r'M$_{infall}$',r'R$_{disrupt}$',r'R$_{merger}$'],
                           quantiles=[0.016,0.5,0.84],
                           title_fmt='.3f',
                           show_titles=True)
    
if(not param.get_bool("fit_r_merger") and not param.get_bool("fit_r_fof")):
    data = np.column_stack((np.log10(mcmc_m_i[slct_half]),mcmc_r_d[slct_half]))
    figure = corner.corner(data,labels=[r'M$_{infall}$',r'R$_{disrupt}$'],
                           quantiles=[0.016,0.5,0.84],
                           title_fmt='.3f',
                           show_titles=True)
    

plt.figure()
max_id = np.max(mcmc_id)
for i in range(0,max_id+1):
    slct = mcmc_id == i
    plt.plot(mcmc_step[slct],mcmc_val[slct],'-',alpha=0.3)
    
plt.ylabel('-Log Likelihood')
plt.xlabel('step number')
plt.yscale('log')
print __file__
dtk.save_figs("figs/"+param.file_name+"/"+__file__+"/",extension='.png')

plt.show()


