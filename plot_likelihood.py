#!/usr/bin/env python2.7
import numpy as np
import matplotlib
import os
#checks if there is a display to use.
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from core_fit2_util import *
from catalog_reader import Catalog,frag_to_real
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize





total_start = time.time()

print "2nd fitting code"
param = dtk.Param(sys.argv[1])

n2lib_loc = param.get_string("n2lib_loc")


fof_loc = param.get_string("fof_loc")
sod_loc = param.get_string("sod_loc")
core_loc = param.get_string("core_loc")
zmr_loc = param.get_string("zmr_loc")
#processed_core_loc = param.get_string('processed_core_loc').replace('${param}',param.file) #output file

fof_cat = Catalog(fof_loc)
sod_cat = Catalog(sod_loc)
core_cat= Catalog(core_loc)
zmr_sdss = np.load(zmr_loc)
zmr_sdss = npzfile_to_dic(zmr_sdss)
print zmr_sdss.keys()

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
core_cat.add_data_name("vx")
core_cat.add_data_name("vy")
core_cat.add_data_name("vz")
core_cat.add_data_name("radius")
core_cat.add_data_name("infall_mass")
core_cat.add_data_name("infall_step")
core_cat.add_data_name("core_tag")
print "reading in files"
core_cat.read_gio()


print "reading in gio files"
fof_cat.read_gio()
sod_cat.read_gio()

print "merging catalogs"
halo_cat = Catalog()
halo_cat.join(fof_cat,sod_cat,join_on='fof_halo_tag')
halo_cat.sort('fof_halo_tag')
core_cat.apply_function('fof_halo_tag',frag_to_real)
clstrs = []
zmr_valid = ZMRIndexValidator(zmr_sdss)
#core_cat.make_kdtree()
avg = []
avg2 = []
for step in steps:
    htags = halo_cat[step]['fof_halo_tag']
    for i in range(0,htags.size):
        if(i%10000==0):
            print i,"/",htags.size
        htag = halo_cat[step]['fof_halo_tag'][i]
        m200 = halo_cat[step]['sod_halo_mass'][i]
        r200 = halo_cat[step]['sod_halo_radius'][i]
        x = halo_cat[step]['fof_halo_center_x'][i]
        y = halo_cat[step]['fof_halo_center_y'][i]
        z = halo_cat[step]['fof_halo_center_z'][i]

        if(m200<1e14):
            continue
        clstr = Cluster(stepz.get_z(step),m200,r200,x,y,z,zmr_valid)
        box = True



        cores_slct = core_cat.cut_box(step,x,y,z,2,2,2,256,256,256,ignore_periodic=False)
        cores_x = core_cat[step]['x'][cores_slct]
        cores_y = core_cat[step]['y'][cores_slct]
        cores_z = core_cat[step]['z'][cores_slct]
        cores_r = core_cat[step]['radius'][cores_slct]
        cores_m = core_cat[step]['infall_mass'][cores_slct]
        #cores2_x,cores2_y,cores2_z,cores2_r = core_cat.cut_box_kdtree(step,x,y,z,2)

        if(box):
            indx = np.arange(cores_x.size)
            numpy.random.shuffle(indx)
            indx = np.argsort(cores_x)
            cores_x = cores_x[indx].astype('float64')
            cores_y = cores_y[indx].astype('float64')
            cores_z = cores_z[indx].astype('float64')
            cores_r = cores_r[indx].astype('float64')
            cores_m = cores_m[indx].astype('float64')
            clstr.set_cores(cores_x,cores_y,cores_z,cores_r,cores_m,ignore_periodic=False)
        else:
            indx = np.arange(cores2_x.size)
            numpy.random.shuffle(indx)
            cores2_x = cores2_x[indx].astype('float64')
            cores2_y = cores2_y[indx].astype('float64')
            cores2_z = cores2_z[indx].astype('float64')
            cores2_r = cores2_r[indx].astype('float64')
            cores2_m = cores2_m[indx].astype('float64')
            clstr.set_cores(cores2_x,cores2_y,cores2_z,cores2_r,cores2_m,ignore_periodic=False)
        
        avg.append(cores_x.size)
        avg2.append(cores_x.size**2)
        clstr.set_n2merg(n2merger)
        clstrs.append(clstr)


print "avg: ",np.average(avg)
print "avg: ",np.average(avg2)
print "done setting the damn things up"



if(True):#make2dplot=='y'):
    #Making the 2d cost plot
    dis_bins = np.logspace(np.log10(0.001),np.log10(1),25)
    dis_avg  = np.sqrt(dis_bins[1:]*dis_bins[:-1])
    #merg_bins = np.logspace(np.log10(0.0001),np.log10(0.5),5)
    merg_bins = np.linspace(0,.15,25)
    merg_avg = np.sqrt(merg_bins[1:]*merg_bins[:-1])
    
    infall_bins = np.linspace(10.0,12,25)
    infall_avg = (infall_bins[1:]+infall_bins[:-1])/2.0
    
    cost_mat_Ngal = np.zeros((dis_avg.size,merg_avg.size))
    cost_mat_gal_den = np.zeros((dis_avg.size,merg_avg.size,infall_avg.size))
    print "starting the param grid search"
    time1 = time.time()
    for i in range(0,dis_avg.size):
        print i,"/",dis_avg.size
        for j in range(0,merg_avg.size):
            for k in range(0,infall_avg.size):
                core_zmr = zmr_from_clusters(dis_avg[i],merg_avg[j],clstrs,zmr_valid,infall_mass=infall_avg[k])
                cost_mat_gal_den[i,j,k]= calc_gal_density_cost2(core_zmr,zmr_sdss,)
        #for k in range(0,merg_avg.size):
        #plt.figure()
        #plt.title("Gal Density Cost, infall=10e%.2f\nmin cost: %.0f"%(infall_avg[k],np.min(cost_mat_gal_den[:,k,:])))
        #plt.pcolor(dis_bins,infall_bins,cost_mat_gal_den[:,k,:].T,norm=LogNorm(vmin=np.min(cost_mat_gal_den),vmax=np.max(cost_mat_gal_den)))
        #plt.xlabel('disruption length [Mpc/h]')
        #plt.ylabel('merger length [Mpc/h]')
        #plt.xscale('log')
        #plt.colorbar()
    likelihood = np.exp(-cost_mat_gal_den)
    likelihood /=np.sum(likelihood)
    dis_collapse = np.sum(likelihood,axis=(1,2))
    plt.figure()
    plt.title('R_disrupt likelihood')
    plt.plot(dis_avg,dis_collapse)
    plt.ylabel('Likelihood')
    plt.xlabel('R_disrupt [h^-1 Mpc]')

    merg_collapse = np.sum(likelihood,axis=(0,2))
    plt.figure()
    plt.title('R_merger likelihood')
    plt.plot(merg_avg,merg_collapse)
    plt.ylabel('Likelihood')
    plt.xlabel('R_merger [h^-1 Mpc]')

    infall_collapse = np.sum(likelihood,axis=(0,1))
    plt.figure()
    plt.title('M_infall likelihood')
    plt.plot(infall_avg,infall_collapse)
    plt.ylabel('Likelihood')
    plt.xlabel('M_infall [h^-1 M_sun]')


    sigmas = (.68,.87)
    dis_infall_collapse = np.sum(likelihood,axis=1)
    lvls = dtk.conf_interval(dis_infall_collapse,sigmas)
    print sigmas
    print lvls
    print np.sum(dis_infall_collapse)
    plt.figure()
    plt.pcolor(dis_bins,infall_bins,dis_infall_collapse.T,norm=LogNorm(),cmap='PuBu')
    #plt.contour(dis_avg,infall_avg,dis_infall_collapse.T,lvls)
    plt.xlabel('R_disrupt [h^-1 Mpc]')
    plt.ylabel('M_infall [h^-1 M_sun]')
    plt.xscale('log')

    dis_merg_collapse = np.sum(likelihood,axis=2)
    lvls = dtk.conf_interval(dis_merg_collapse,sigmas)
    plt.figure()
    plt.pcolor(dis_bins,merg_bins,dis_merg_collapse.T,norm=LogNorm(),cmap='PuBu')#,norm=LogNorm(),
    plt.contour(dis_avg,merg_avg,dis_merg_collapse.T,lvls)

    plt.xlabel('R_disrupt [h^-1 Mpc]')
    #plt.xlabel('M_infall [h^-1 M_sun]')
    plt.ylabel('R_merger [h^-1 Mpc]')
    plt.xscale('log')

    merg_infall_collapse = np.sum(likelihood,axis=0)
    lvls = dtk.conf_interval(merg_infall_collapse,sigmas)
    plt.figure()
    plt.pcolor(merg_bins,infall_bins,merg_infall_collapse.T,norm=LogNorm(),cmap='PuBu')#,norm=LogNorm(),
    plt.contour(merg_avg,infall_avg,merg_infall_collapse.T,lvls)
    #plt.xlabel('R_disrupt [h^-1 Mpc]')
    plt.ylabel('M_infall [h^-1 M_sun]')
    plt.xlabel('R_merger [h^-1 Mpc]')
    #plt.xscale('log')

    
total_end = time.time()
    
print "total time: ", total_end - total_start
dtk.save_figs("figs/"+param.file+"/plot_likelihood.py/",extension='.png')

plt.show()

