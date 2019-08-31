#!/usr/bin/env python2.7

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from zmr import ZMR
import dtk
import sys 

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size': 12,'weight':'bold'})
rc('text',usetex=True)


param_file1 = sys.argv[1]
param_file2 = sys.argv[2]
if(len(sys.argv) == 5):
    param_label1 = sys.argv[3]
    param_label2 = sys.argv[4]
else:
    param_label1 = sys.argv[1]
    param_label2 = sys.argv[2]
zmr_cores1 = ZMR(file_loc="output/"+param_file1+"/zmr_lkhd_cores.param")
zmr_cores2 = ZMR(file_loc="output/"+param_file2+"/zmr_lkhd_cores.param")
#zmr_sdss_npz = np.load("/home/dkorytov/phys/Ngal_sdss/data/normal_mask4/result/type1_weight1_mag1_clr1_result.npz")
#print zmr_sdss_npz.keys()

#area = np.pi*(zmr_sdss.r_bins[1:]**2-zmr_sdss.r_bins[:-1]**2)
#for zi in range(0,zmr_sdss.z_size):
#    for mi in range(0,zmr_sdss.z_size):
#        zmr_cores.zmr_gal_density[zi,mi,:]=zmr_cores.zmr_gal_counts[zi,mi,:]/area/zmr_cores.zm_counts[zi,mi]
        

r_avg = (zmr_cores1.r_bins[:-1]+zmr_cores1.r_bins[1:])/2.0

colors = ['b','r','g','m','c','y']
c_i = 0


for zi in range(0,zmr_cores2.z_bins.size-1):
    if(np.sum(zmr_cores2.zm_counts[zi,:])==0 or np.sum(zmr_cores1.zm_counts[zi,:])==0):
        continue #Both don't have data here
    plt.figure()
    plt.title("Galaxy Density %.2f $<$ z $<$ %.2f"%(zmr_cores2.z_bins[zi],zmr_cores2.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_cores2.m_bins.size-1):
        if(zmr_cores2.zm_counts[zi,mi] > 0 and zmr_cores1.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gd = zmr_cores2.zmr_gal_density[zi,mi]
            err = zmr_cores2.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            plt.plot(r_avg,zmr_gd,'s--',color=c,lw=2,mfc='none',mec=c,mew=1.5,ms=8)
            #plt.scatter(r_avg,zmr_gd,marker='s',facecolor='none',edgecolor=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            zmr_gd = zmr_cores1.zmr_gal_density[zi,mi]
            err = zmr_cores1.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            zmr_gd = np.clip(zmr_gd,1e-1,np.max(zmr_gd))
            plt.plot(r_avg,zmr_gd,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            plt.plot([],[],label=r'%.2e$<$M200$<$%.2e'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1]),color=c)
    plt.plot([],[],'ks--', lw=2, label=param_label2, mfc='none', mew=1.5, ms=8)
    plt.plot([],[],'ko-', label=param_label1)
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$r/R_{200}$')
    plt.ylabel(r'$\Sigma_{gal}$[$gal/R_{200}^{2}$]')
    plt.yscale('log')

for mi in range(0,zmr_cores2.m_bins.size-1):
    if(np.sum(zmr_cores2.zm_counts[:,mi])==0 or np.sum(zmr_cores1.zm_counts[:,mi])==0):
        continue #Both don't have data here
    plt.figure()
    #plt.title(r"Galaxy Density %.2f $<$ z $<$ %.2f"%(zmr_cores2.z_bins[zi],zmr_cores2.z_bins[zi+1]))
    plt.title(r'%.2e$<$M200$<$%.2e'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1]))
    c_i=0
    for zi in range(0,zmr_cores2.z_bins.size-1):
        if(zmr_cores2.zm_counts[zi,mi] > 0 and zmr_cores1.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gd = zmr_cores2.zmr_gal_density[zi,mi]
            err = zmr_cores2.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            plt.plot(r_avg,zmr_gd,'s--',color=c,lw=2,mfc='none',mec=c,mew=1.5,ms=8)
            #plt.scatter(r_avg,zmr_gd,marker='s',facecolor='none',edgecolor=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            zmr_gd = zmr_cores1.zmr_gal_density[zi,mi]
            err = zmr_cores1.zmr_gal_density_err[zi,mi]
            min_err = np.clip(zmr_gd-err,1e-1,np.max(zmr_gd-err))
            zmr_gd = np.clip(zmr_gd,1e-1,np.max(zmr_gd))
            plt.plot(r_avg,zmr_gd,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_gd+err,color=c,alpha=0.3)
            plt.plot
            #plt.plot([],[],label=r'%.2e$<$M200$<$%.2e'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1]),color=c)
            plt.plot
            plt.plot([],[],label=r"%.2f $<$ z $<$ %.2f"%(zmr_cores2.z_bins[zi],zmr_cores2.z_bins[zi+1]))
    counts1 = int(np.sum(zmr_cores1.zm_counts[:,mi]))
    counts2 = int(np.sum(zmr_cores2.zm_counts[:,mi]))
    plt.plot([],[],'ks--',lw=2,label=param_label2.replace("_","\_")+"[{:}]".format(counts2),mfc='none',mew=1.5,ms=8)
    plt.plot([],[],'ko-',label=param_label1.replace("_","\_")+"[{:}]".format(counts1))
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$r/R_{200}$')
    plt.ylabel(r'$\Sigma_{gal}$[$gal/R_{200}^{2}$]')
    plt.yscale('log')

for zi in range(0,zmr_cores2.z_bins.size-1):
    if(np.sum(zmr_cores2.zm_counts[zi,:])==0 or np.sum(zmr_cores1.zm_counts[zi,:])==0):
        continue #Both don't have data here
    print(zmr_cores1.zm_counts[zi,:])
    print(zmr_cores2.zm_counts[zi,:])
    plt.figure()
    plt.title("Galaxy dn/dr \n%.2f$<$z$<$%.2f"%(zmr_cores2.z_bins[zi],zmr_cores2.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_cores2.m_bins.size-1):
        if(zmr_cores2.zm_counts[zi,mi] > 0 and zmr_cores1.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_dgdr = zmr_cores2.zmr_dgal_dr[zi,mi]
            err = zmr_cores2.zmr_dgal_dr_err[zi,mi]
            min_err = np.clip(zmr_dgdr-err,1e-1,np.max(zmr_dgdr-err))
            plt.plot(r_avg,zmr_dgdr,'x--',color=c)
            plt.fill_between(r_avg,min_err,zmr_dgdr+err,color=c,alpha=0.5)
            zmr_dgdr = zmr_cores1.zmr_dgal_dr[zi,mi]
            err = zmr_cores1.zmr_dgal_dr_err[zi,mi]
            min_err = np.clip(zmr_dgdr-err,1e-1,np.max(zmr_dgdr-err))
            zmr_dgdr = np.clip(zmr_dgdr,1e-1,np.max(zmr_dgdr))
            plt.plot(r_avg,zmr_dgdr,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_dgdr+err,color=c,alpha=0.5)
            plt.plot([],[],label='%.2e$<$M200$<$%.2e'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1]),color=c)
    plt.plot([],[],'kx--',label=param_label2.replace("_","\_")+"[{}]".format(counts2))
    plt.plot([],[],'ko-',label=param_label1.replace("_","\_")+"[{}]".format(counts1))
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'$dgal/dr [1/Mpc h^{-1}]$')
    plt.yscale('log')

for zi in range(0,zmr_cores2.z_bins.size-1):
    if(np.sum(zmr_cores2.zm_counts[zi,:])==0 or np.sum(zmr_cores1.zm_counts[zi,:])==0):
        continue #Both don't have data here
    print(zmr_cores1.zm_counts[zi,:])
    print(zmr_cores2.zm_counts[zi,:])
    plt.figure()
    plt.title("Galaxy Accumulated Count \n%.2f$<$z$<$%.2f"%(zmr_cores2.z_bins[zi],zmr_cores2.z_bins[zi+1]))
    c_i=0
    for mi in range(0,zmr_cores2.m_bins.size-1):
        if(zmr_cores2.zm_counts[zi,mi] > 0 and zmr_cores1.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gacm = zmr_cores2.zmr_gal_accum[zi,mi]
            err = zmr_cores2.zmr_gal_accum_err[zi,mi]
            min_err = np.clip(zmr_gacm-err,1e-1,np.max(zmr_gacm-err))
            plt.plot(r_avg,zmr_gacm,'s--',lw=3,color=c)
            plt.fill_between(r_avg,min_err,zmr_gacm+err,color=c,alpha=0.5)
            zmr_gacm = zmr_cores1.zmr_gal_accum[zi,mi]
            err = zmr_cores1.zmr_gal_accum_err[zi,mi]
            min_err = np.clip(zmr_gacm-err,1e-1,np.max(zmr_gacm-err))
            zmr_gacm = np.clip(zmr_gacm,1e-1,np.max(zmr_gacm))
            plt.plot(r_avg,zmr_gacm,'o-',color=c)
            plt.fill_between(r_avg,min_err,zmr_gacm+err,color=c,alpha=0.5)
            plt.plot([],[],label='%.2e$<$M200$<$%.2e'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1]),color=c)
    plt.plot([],[],'ks--',label=param_label2.replace("_","\_")+"[{}]".format(counts2))
    plt.plot([],[],'ko-',label=param_label1.replace("_","\_")+"[{}]".format(counts1))
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'$Ngal(r)$')
    plt.yscale('log')


for zi in range(0,zmr_cores2.z_bins.size-1):
    if(np.sum(zmr_cores2.zm_counts[zi,:])==0 or np.sum(zmr_cores1.zm_counts[zi,:])==0):
        continue #at least doesn't have data here
    plt.figure()
    c_i =0;
    total_sum = 0
    for mi in range(0,zmr_cores2.m_bins.size-1):
        if(zmr_cores2.zm_counts[zi,mi] > 0 and zmr_cores1.zm_counts[zi,mi] > 0):
            c = colors[c_i%len(colors)]
            c_i +=1
            zmr_gd_sdss = zmr_cores2.zmr_gal_density[zi,mi]
            zmr_err_sdss = zmr_cores2.zmr_gal_density_err[zi,mi]
            zmr_gd_core = zmr_cores1.zmr_gal_density[zi,mi]
            zmr_err_core = zmr_cores1.zmr_gal_density_err[zi,mi]
            diff = (zmr_gd_sdss - zmr_gd_core)**2
            err  = zmr_err_sdss**2 + zmr_err_core**2
            res = diff/err
            total = np.sum(res)
            total_sum += np.sum(res)
            plt.plot(r_avg,res,'o-',color=c,label='%.2e$<$M200$<$%.2e:%f'%(zmr_cores2.m_bins[mi],zmr_cores2.m_bins[mi+1],total))
            if(c_i ==9):
                for i in range(0,10):
                    print( "diff:")
                    print( zmr_gd_sdss[i],"-",zmr_gd_core[i],"=", zmr_gd_sdss[i]-zmr_gd_core[i])
                    print( "->",diff[i])
                    print( "err:")
                    print( zmr_err_sdss[i],"+",zmr_err_core[i],"=",err[i])
                    print( "res:")
                    print( res[i])
    plt.title("Source of Error [%f]"%total_sum)
    plt.grid()
    plt.legend(loc='best',framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'\text{Error}')
    plt.yscale('log')
#dtk.set_fig_path("figs/zmrs/")
dtk.save_figs("figs/"+param_file1+"/"+__file__+"/"+param_file2+"/",extension=".png")
plt.show()
