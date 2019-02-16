#!/usr/bin/env python2.7

import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from zmr import ZMR
import dtk
import sys 

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size': 12,'weight':'bold'})
rc('text',usetex=True)

param_file = sys.argv[1]

zmr_sdss = ZMR(file_loc="output/"+param_file+"/zmr_sdss.param")
if dtk.file_exists("output/"+param_file+"/zmr_lkhd_cores.param"):
    zmr_cores = ZMR(file_loc="output/"+param_file+"/zmr_lkhd_cores.param")
    print "likelihood zmrs"
else:
    zmr_cores = ZMR(file_loc="output/"+param_file+"/zmr_cores.param")
    print "fit zmrs"
#zmr_sdss_npz = np.load("/home/dkorytov/phys/Ngal_sdss/data/normal_mask4/result/type1_weight1_mag1_clr1_result.npz")
#print zmr_sdss_npz.keys()

#area = np.pi*(zmr_sdss.r_bins[1:]**2-zmr_sdss.r_bins[:-1]**2)
#for zi in range(0,zmr_sdss.z_size):
#    for mi in range(0,zmr_sdss.z_size):
#        zmr_cores.zmr_gal_density[zi,mi,:]=zmr_cores.zmr_gal_counts[zi,mi,:]/area/zmr_cores.zm_counts[zi,mi]
        

r_avg = (zmr_sdss.r_bins[:-1]+zmr_sdss.r_bins[1:])/2.0

colors = ['b','r','g','m','c','y']
c_i = 0


#Plot all overlapping mass bins
for zi in range(0,zmr_sdss.z_bins.size-1):
    if(np.sum(zmr_sdss.zm_counts[zi,:])==0 or np.sum(zmr_cores.zm_counts[zi,:])==0):
        continue #Both don't have data here
    plt.figure()
    plt.title(r"Galaxy Density %.2f $<$ z $<$ %.2f"%(zmr_sdss.z_bins[zi],zmr_sdss.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_sdss.m_bins.size-1):
        if(zmr_sdss.zm_counts[zi,mi] > 0 and zmr_cores.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gd = zmr_sdss.zmr_gal_density[zi,mi]
            err = zmr_sdss.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            plt.plot(r_avg,zmr_gd,'s--',color=c,lw=2,mfc='none',mec=c,mew=1.5,ms=8)
            #plt.scatter(r_avg,zmr_gd,marker='s',facecolor='none',edgecolor=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            zmr_gd = zmr_cores.zmr_gal_density[zi,mi]
            err = zmr_cores.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            zmr_gd = np.clip(zmr_gd,1e-1,np.max(zmr_gd))
            plt.plot(r_avg,zmr_gd,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            plt.plot([],[],label=r'%.2e$<$M200$<$%.2e'%(zmr_sdss.m_bins[mi],zmr_sdss.m_bins[mi+1]),color=c)
    plt.plot([],[],'ks--',lw=2,label=r'RedMapper Gal. Profile',mfc='none',mew=1.5,ms=8)
    plt.plot([],[],'ko-',label=r'Core Profile')
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$r/R_{200}$')
    plt.ylabel(r'$\Sigma_{gal}$[$gal/R_{200}^{2}$]')
    plt.yscale('log')
c_i=0

#plot individual mass bins
gs=gridspec.GridSpec(5,1)
#gs.update(hspace=0.05)
for mi in range(0,zmr_sdss.m_bins.size-1):
    if(np.sum(zmr_sdss.zm_counts[:,mi])==0 or np.sum(zmr_cores.zm_counts[:,mi])==0):
        continue #Both don't have data here
    plt.figure()
    #plt.title(r"Galaxy Density %.2f $<$ z $<$ %.2f"%(zmr_sdss.z_bins[zi],zmr_sdss.z_bins[zi+1]))
    plt.title(r'%.2e$<$M200$<$%.2e'%(zmr_sdss.m_bins[mi],zmr_sdss.m_bins[mi+1]))
    c = colors[c_i%len(colors)]
    c_i +=1
    for zi in range(0,zmr_sdss.z_bins.size-1):
        if(zmr_sdss.zm_counts[zi,mi] > 0 and zmr_cores.zm_counts[zi,mi] > 0):
            ax1 = plt.subplot(gs[:4,:])
            zmr_plt_sdss = zmr_sdss.zmr_gal_density[zi,mi]
            err = zmr_sdss.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_plt_sdss-err,1e-1,np.max(zmr_plt_sdss-err))
            ax1.plot(r_avg,zmr_plt_sdss,'s--',color=c,lw=2,mfc='none',mec=c,mew=1.5,ms=8)
            #plt.scatter(r_avg,zmr_plt_sdss,marker='s',facecolor='none',edgecolor=c)
            ax1.fill_between(r_avg,min_err,zmr_plt_sdss+err,color=c,alpha=0.3)
            zmr_plt_cores = zmr_cores.zmr_gal_density[zi,mi]
            err = zmr_cores.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_plt_cores-err,1e-1,np.max(zmr_plt_cores-err))
            zmr_plt_cores = np.clip(zmr_plt_cores,1e-1,np.max(zmr_plt_cores))
            ax1.plot(r_avg,zmr_plt_cores,'o-',color=c)
            ax1.fill_between(r_avg,min_err,zmr_plt_cores+err,color=c,alpha=0.3)
            plt.yscale('log')
    
            #plt.plot([],[],label=r'%.2e$<$M200$<$%.2e'%(zmr_sdss.m_bins[mi],zmr_sdss.m_bins[mi+1]),color=c)

            plt.plot([],[],label=r"%.2f $<$ z $<$ %.2f"%(zmr_sdss.z_bins[zi],zmr_sdss.z_bins[zi+1]))
            ax2 = plt.subplot(gs[4,:])
            yerr = np.sqrt( (zmr_sdss.zmr_gal_density_err[zi,mi]/zmr_plt_sdss)**2 + (zmr_cores.zmr_gal_density_err[zi,mi]/zmr_plt_cores)**2)
            ax2.errorbar(r_avg, zmr_plt_cores/zmr_plt_sdss,  yerr = yerr, fmt = 'o', color=c)
            ax2.set_ylim((0.0,2.0))
            ax2.axhline(1.0, color='k', ls='--')
            ax2.set_ylabel('$\Sigma_{cores}/\Sigma_{sdss}$')
    ax1.plot([],[],'ks--',lw=2,label=r'RedMapper Gal. Profile',mfc='none',mew=1.5,ms=8)
    ax1.plot([],[],'ko-',label=r'Core Profile')
    ax1.get_xaxis().set_visible(False)
    ax1.grid()
    ax1.legend(loc='best',framealpha=0.5)
    ax2.set_xlabel(r'$r/R_{200}$')
    ax1.set_ylabel(r'$\Sigma$[$galaxies/R_{200}^{2}$]')


for zi in range(0,zmr_sdss.z_bins.size-1):
    if(np.sum(zmr_sdss.zm_counts[zi,:])==0 or np.sum(zmr_cores.zm_counts[zi,:])==0):
        continue #Both don't have data here
    print zmr_cores.zm_counts[zi,:]
    print zmr_sdss.zm_counts[zi,:]
    plt.figure()
    plt.title("Galaxy dn/dr \n%.2f$<$z$<$%.2f"%(zmr_sdss.z_bins[zi],zmr_sdss.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_sdss.m_bins.size-1):
        if(zmr_sdss.zm_counts[zi,mi] > 0 and zmr_cores.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_dgdr = zmr_sdss.zmr_dgal_dr[zi,mi]
            err = zmr_sdss.zmr_dgal_dr_err[zi,mi]
            min_err = np.clip(zmr_dgdr-err,1e-1,np.max(zmr_dgdr-err))
            plt.plot(r_avg,zmr_dgdr,'x--',color=c)
            plt.fill_between(r_avg,min_err,zmr_dgdr+err,color=c,alpha=0.5)
            zmr_dgdr = zmr_cores.zmr_dgal_dr[zi,mi]
            err = zmr_cores.zmr_dgal_dr_err[zi,mi]
            min_err = np.clip(zmr_dgdr-err,1e-1,np.max(zmr_dgdr-err))
            zmr_dgdr = np.clip(zmr_dgdr,1e-1,np.max(zmr_dgdr))
            plt.plot(r_avg,zmr_dgdr,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_dgdr+err,color=c,alpha=0.5)
            plt.plot([],[],label='%.2e$<$M200$<$%.2e'%(zmr_sdss.m_bins[mi],zmr_sdss.m_bins[mi+1]),color=c)
    plt.plot([],[],'kx--',label='RedMapper Gal. Profile')
    plt.plot([],[],'ko-',label='Core Profile')
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'$dgal/dr [1/Mpc h^{-1}]$')
    plt.yscale('log')

for zi in range(0,zmr_sdss.z_bins.size-1):
    if(np.sum(zmr_sdss.zm_counts[zi,:])==0 or np.sum(zmr_cores.zm_counts[zi,:])==0):
        continue #Both don't have data here
    print zmr_cores.zm_counts[zi,:]
    print zmr_sdss.zm_counts[zi,:]
    plt.figure()
    plt.title("Galaxy Accumulated Count \n%.2f$<$z$<$%.2f"%(zmr_sdss.z_bins[zi],zmr_sdss.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_sdss.m_bins.size-1):
        if(zmr_sdss.zm_counts[zi,mi] > 0 and zmr_cores.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gacm = zmr_sdss.zmr_gal_accum[zi,mi]
            err = zmr_sdss.zmr_gal_accum_err[zi,mi]
            min_err = np.clip(zmr_gacm-err,1e-1,np.max(zmr_gacm-err))
            plt.plot(r_avg,zmr_gacm,'s--',lw=3,color=c)
            plt.fill_between(r_avg,min_err,zmr_gacm+err,color=c,alpha=0.5)
            zmr_gacm = zmr_cores.zmr_gal_accum[zi,mi]
            err = zmr_cores.zmr_gal_accum_err[zi,mi]
            min_err = np.clip(zmr_gacm-err,1e-1,np.max(zmr_gacm-err))
            zmr_gacm = np.clip(zmr_gacm,1e-1,np.max(zmr_gacm))
            plt.plot(r_avg,zmr_gacm,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_gacm+err,color=c,alpha=0.5)
            plt.plot([],[],label='%.2e$<$M200$<$%.2e'%(zmr_sdss.m_bins[mi],zmr_sdss.m_bins[mi+1]),color=c)
    plt.plot([],[],'ks--',label='RedMapper Gal. Profile')
    plt.plot([],[],'ko-',label='Core Profile')
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'$Ngal(r)$')
    plt.yscale('log')


for zi in range(0,zmr_sdss.z_bins.size-1):
    if(np.sum(zmr_sdss.zm_counts[zi,:])==0 or np.sum(zmr_cores.zm_counts[zi,:])==0):
        continue #at least doesn't have data here
    plt.figure()
    c_i =0;
    total_sum = 0
    dof = 0
    for mi in range(0,zmr_sdss.m_bins.size-1):
        if(zmr_sdss.zm_counts[zi,mi] > 0 and zmr_cores.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gd_sdss = zmr_sdss.zmr_gal_density[zi,mi]
            zmr_err_sdss = zmr_sdss.zmr_gal_density_err[zi,mi]
            zmr_gd_core = zmr_cores.zmr_gal_density[zi,mi]
            zmr_err_core = zmr_cores.zmr_gal_density_err[zi,mi]
            diff = (zmr_gd_sdss - zmr_gd_core)**2
            err  = zmr_err_sdss**2 + zmr_err_core**2
            res = diff/err
            total = np.sum(res)
            total_sum += np.sum(res)
            dof += np.size(res)
            plt.plot(r_avg, res, '-',color=c,label='%.2f$<$M200$<$%.2f err:%.1f'%(np.log10(zmr_sdss.m_bins[mi]),np.log10(zmr_sdss.m_bins[mi+1]),total))
            plt.fill_between(r_avg, 0, res, alpha=0.3, color=c)
            if(c_i ==9):
                for i in range(0,10):
                    print "diff:"
                    print zmr_gd_sdss[i],"-",zmr_gd_core[i],"=", zmr_gd_sdss[i]-zmr_gd_core[i]
                    print "->",diff[i]
                    print "err:"
                    print zmr_err_sdss[i],"+",zmr_err_core[i],"=",err[i]
                    print "res:"
                    print res[i]
    ylim = plt.ylim()
    plt.ylim([0, ylim[1]])
    plt.title("Source of Error [%f] $\chi^{2}$=%f"%(total_sum,(total_sum/2.0)/(dof-2)))
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'\text{Error}')
    # plt.yscale('log')
dtk.set_fig_path("figs/zmrs/")
dtk.save_figs("figs/"+param_file+"/"+__file__+"/",extension=".png")
plt.show()
