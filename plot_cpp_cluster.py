#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import sys
import dtk

param_file_name = sys.argv[1]
i=0
dtk.set_fig_path("figs/"+param_file_name+"/"+__file__+"/")
while(True):

    param = dtk.Param("output/"+param_file_name+"/"+__file__+"/lgrid.param")
    # if(sys.argv[1] == 'cpp'):
    #     param = dtk.Param("tmp/clstr%d.param"%i)
    # elif(sys.argv[1] == 'py'):
    #     param = dtk.Param("tmp/clstr%d_py.param"%i)
    # else:
    #     print "no \"cpp\"/\"py\" argument afterwards. Quitting"
    halo_x = param.get_float("halo_x")
    halo_y = param.get_float("halo_y")
    halo_m = param.get_float("halo_mass")
    halo_r = param.get_float("halo_radius")

    core_x = np.array(param.get_float_list("core_x"))
    core_y = np.array(param.get_float_list("core_y"))
    core_r = np.array(param.get_float_list("core_r"))
    core_m = np.array(param.get_float_list("core_m"))
    print np.min(core_r),"->",np.max(core_r)
    cprtcl_x = np.array(param.get_float_list("cprtcl_x"))
    cprtcl_y = np.array(param.get_float_list("cprtcl_y"))

    dis_cprtcl_x = np.array(param.get_float_list("dis_cp_x"))    
    dis_cprtcl_y = np.array(param.get_float_list("dis_cp_y"))    

    gal_x = np.array(param.get_float_list("gal_x"))
    gal_y = np.array(param.get_float_list("gal_y"))
    gal_w = np.array(param.get_float_list("gal_w"))
    gal_type = np.array(param.get_float_list("gal_type"))
    
    r_bins = np.array(param.get_float_list("r_bins"))
    r_cnt = np.array(param.get_float_list("r_cnt"))
    ngal = param.get_float("Ngal")
    m_infall = param.get_float("m_infall")
    r_disrupt = param.get_float("r_disrupt")
    r_merger = param.get_float("r_merger")
#    if(halo_m < 7e14):
#        i +=1
#        continue
#    if( halo_m > 7.74e14 or halo_m < 4.64e14):
#        i+=1
#        continue
    gal_slct_r  = (core_r < r_disrupt) 
    gal_slct_m  = (core_m > m_infall)
    gal_slct    = gal_slct_m & gal_slct_r
    #print sys.argv[1],m_infall,r_disrupt,r_merger,"from ",i," --> ",np.sum(gal_slct),"/",np.sum(gal_slct_m),'/',np.sum(gal_slct_r),'/',core_r.size
    indx = np.argmax(core_m)
    #print sys.argv[1],core_m[indx],core_r[indx]
    #plt.figure()
    #plt.title(sys.argv[1])
    #plt.hist(core_r)

    #plt.figure(figsize=(12,6))
    plt.figure(figsize=(24,12))
    ax0 = plt.subplot(121)
    plt.title("SOD mass: %.2e Msun/h"%(halo_m))
    plt.plot(cprtcl_x,cprtcl_y,'g+',label='all core prtcls')
    plt.plot(core_x,core_y,'s',mec='r',mfc='none',mew=1,label='all cores[%d]'%core_x.size)
    plt.plot(halo_x,halo_y,'cs',label='halo center')
    plt.plot([],[],'k',label='SOD radius'%halo_r)
    plt.xlabel('x [h^-1 Mpc]')
    plt.ylabel('y [h^-1 Mpc]')
    plt.axis('equal')
    plt.grid()
    plt.legend(loc='best',framealpha=.5)
    for j in range(0,r_bins.size):
        circle0 = plt.Circle((halo_x,halo_y),halo_r*r_bins[j],facecolor='none',edgecolor='k')
        plt.gca().add_artist(circle0)
    slct_gal = gal_type == 0
    slct_bcg = gal_type == 1
    ax1 = plt.subplot(122,sharex=ax0,sharey=ax0)
    for j in range(0,r_bins.size):
        circle1 = plt.Circle((halo_x,halo_y),halo_r*r_bins[j],facecolor='none',edgecolor='k')
        ax1.add_artist(circle1)
    plt.title("SOD radius:%.2f"%(halo_r))
    g_alpha = clr.colorConverter.to_rgba('green',alpha=.0)
    m_alpha = clr.colorConverter.to_rgba('magenta',alpha=1.0)
    plt.scatter(dis_cprtcl_x,dis_cprtcl_y,marker='x',lw=0.01,c='g',label='disrupted core prtcls')

 

    plt.scatter(gal_x[slct_gal],gal_y[slct_gal],facecolors='none',edgecolors='r',s=gal_w[slct_gal]*20,
                label='compact cores[%d/%d/%d]'%(np.sum(slct_gal),np.sum(gal_w[slct_gal]),core_x.size))
    plt.scatter(gal_x[slct_bcg],gal_y[slct_bcg],marker='o',color=m_alpha,s=30,
                label='disrupted fof groups[%d]'%np.sum(slct_bcg))
    plt.plot(halo_x,halo_y,'s',label='halo center',mec='c',mfc='none',mew=1.5)
    plt.plot([],[],'k',label='SOD radius'%halo_r)
    plt.xlabel('x [h^-1 Mpc]')
    plt.ylabel('y [h^-1 Mpc]')
    plt.axis('equal')
    plt.grid()

    if(np.sum(slct_gal )+np.sum(slct_bcg) > 0):
        plt.legend(loc='best',framealpha=.5)
    plt.tight_layout()
    #plt.show()
    dtk.save_figs(extension='.png')
    plt.close()
    print i
    i+=1
    
