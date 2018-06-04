#!/usr/bin/env python2.7

import numpy as np
from catalog_reader import Catalog,frag_to_real
from core_fit2_util import *
from core_fit2_mcmc import *
import dtk
import sys
import time
import numpy.random
from matplotlib.colors import LogNorm
from scipy.optimize import minimize
from halotools.mock_observables import return_xyz_formatted_array,tpcf

#from astropy.halotools import 

#plt.rc('font', family='san-serif')
#plt.rc('text', usetex=True)


print "running mcmc"
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
start_total = time.time()
core_cat.read_gio()


dis_len = 0.3
infall_mass = 11.6
merg_len =0.08

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

#Add each cluster above 1e14 mass to the cluster list
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
        import pickle
        clstr.set_n2merg(n2merger)
        clstrs.append(clstr)






print "avg: ",np.average(avg)
print "avg: ",np.average(avg2)
print "done setting the damn things up"

zmr_core = zmr_valid.make_empty_zmr()

def zmr_to_min(param):
    dis_len = param[0]
    merg_len = param[1]
    zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid)
    cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
    return cost

#taking infall mass as a variable
def zmr_to_min2(param):
    dis_len = 10**param[0]
    merg_len = param[1]
    infall_m = param[2]
    zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid,infall_m)
    cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
    return cost

end_total = time.time()
print "total process: ",end_total-start_total
start_total = time.time()
x0 = [.3,.081,11.60]
res = minimize(zmr_to_min2,x0,method='nelder-mead',
               options={'xtol':1e-5,'disp':True}) 

print "result: "
print res
print res.x
res.x[0] = np.log10(res.x[0])
dis_len = np.log10(res.x[0])
merg_len = res.x[1]
infall_mass = res.x[2]
    
import corner
print "Running MCMC now..."
chain_len = 1000
nwalker  = 100
start_total = time.time()
sampler = run_mcmc(res.x,nwalker,chain_len,zmr_sdss,zmr_valid,clstrs)

samples = sampler.chain[:,chain_len/2:,:].reshape((-1,3))
np.save("mcmc.npy",samples)
fig = corner.corner(samples, labels=['dis_len','merg_len','infall_m'])


end_total = time.time()
zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid,infall_mass)
print "fit time: ",end_total-start_total


cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
print "this is the cost: ", cost


#Making the 2d cost plot
dis_bins = np.logspace(np.log10(0.01),np.log10(10),2)
dis_avg  = np.sqrt(dis_bins[1:]*dis_bins[:-1])
#merg_bins = np.logspace(np.log10(0.0001),np.log10(0.5),5)
merg_bins = np.linspace(0,.1,2)
merg_avg = np.sqrt(merg_bins[1:]*merg_bins[:-1])

infall_bins = np.linspace(10.0,13,2)
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
for k in range(0,merg_avg.size):
    plt.figure()
    plt.title("Gal Density Cost, infall=10e%.2f\nmin cost: %.0f"%(infall_avg[k],np.min(cost_mat_gal_den[:,k,:])))
    plt.pcolor(dis_bins,infall_bins,cost_mat_gal_den[:,k,:].T,norm=LogNorm(vmin=np.min(cost_mat_gal_den),vmax=np.max(cost_mat_gal_den)))
    plt.xlabel('disruption length [Mpc/h]')
    plt.ylabel('merger length [Mpc/h]')
    plt.xscale('log')
    plt.colorbar()



minimum = 0.1
rad_bin_area = zmr_valid.get_rad_bin_area()
rad_bin_avg = zmr_valid.r_bin_avg
rad_bin = zmr_valid.r_bins

zmr_gal_cnt = zmr_core['zmr_gal_cnt']/rad_bin_area+minimum
zmr_gal_cnt_err = np.sqrt(zmr_core['zmr_gal_cnt'])/rad_bin_area

for i in range(0,zmr_core['z_bins'].size-1):
    for j in range(0,zmr_core['mass_bins'].size-1):
        for k in range(0,zmr_core['rad_bins'].size-1):
            if(zmr_gal_cnt_err[i,j,k] == 0):
                zmr_gal_cnt_err[i,j,k]=1.0

zmr_count_to_density(zmr_core,0)

colors = ['b','r','g','m','c','y']
c_i = 0
for i in range(0,zmr_core['z_bins'].size-1):
    if(np.sum(zmr_core['zm_counts'][i])==0):
        continue
    plt.figure()
    plt.title('Galaxy Clusters Profiles\n%.2f < z < %.2f'%(zmr_sdss['z_bins'][i],zmr_sdss['z_bins'][i+1]))
    c_i=0
    for j in range(0,zmr_core['mass_bins'].size-1):
        if(zmr_core['zm_counts'][i,j]>0 and zmr_sdss['zmr_cnt'][i,j,0]):
            rad_den = zmr_core['zmr_gal_density'][i,j]
            rad_den_err = zmr_core['zmr_gal_density_err'][i,j]
            c = colors[c_i%len(colors)]
            c_i +=1
            plt.plot(rad_bin_avg,zmr_core['zmr_gal_density'][i,j],'o-', color=c)
            plt.fill_between(rad_bin_avg,np.clip(rad_den-rad_den_err,a_min=minimum,a_max=max(rad_den)),rad_den+rad_den_err,alpha='0.5',color=c)
            rad_den = zmr_sdss['rad_prof'][i,j]
            rad_den_err = zmr_sdss['rad_prof_err'][i,j]
            plt.plot(rad_bin_avg,rad_den,'x--',color=c)
            plt.fill_between(rad_bin_avg,np.clip(rad_den-rad_den_err,a_min=minimum,a_max=max(rad_den)),rad_den+rad_den_err,color=c,alpha='0.5')
            plt.plot([],[],color=c,label='%.1e<M200<%.1e'%(zmr_sdss['mass_bins'][j],zmr_sdss['mass_bins'][j+1]))
    plt.plot([],[],'kx--',label='RedMapper Gal. Profile')
    plt.plot([],[],'ko-',label='Core Gal. Profile')
    plt.grid()
    plt.legend(loc='best',fancybox=True,framealpha=0.5)
    plt.xlabel(r'$R/R_{200}$')
    plt.ylabel(r'$\Sigma_{gal}$[$gal/R_{200}^{2}$]')
    plt.yscale('log')
    plt.gca().set_ylim(bottom=1)





dtk.save_figs("figs/"+param.file+"/core_fit2.py/")

plt.show()
