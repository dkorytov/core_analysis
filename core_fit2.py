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


print("2nd fitting code")
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
print(zmr_sdss.keys())

step1 = param.get_int("step")
steps = param.get_int_list("steps")

z_in = param.get_float("z_in")
z_out = param.get_float("z_out")
num_steps = param.get_int("num_steps")
cluster_radial_volume = param.get_float("cluster_radial_volume")
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
print("reading in files")
core_cat.read_gio()


dis_len = 0.3
infall_mass = 11.6
merg_len =0.08
#print "making processed core catalog"
for step in steps:
    continue;
    print("nans x? :",np.sum(np.isnan(core_cat[step]['x'])))
    print("nans radius? :",np.sum(np.isnan(core_cat[step]['radius'])))
    print("nans vx? :",np.sum(np.isnan(core_cat[step]['vx'])))
    print("nans x & radius: ",np.sum(np.isnan(core_cat[step]['radius']) & np.isnan(core_cat[step]['x'])))
    print(np.nanmax(core_cat[step]['radius']))
    print("starting n2 writeout")
    gal_r = core_cat[step]['radius']
    gal_m = core_cat[step]['infall_mass']
    slct = (gal_r < dis_len )& (gal_m > 10**infall_mass)

#    n2merger.n2merger_writeout(step, merg_len,
#                               "/home/dkorytov/proj/core_tracking/core_fit/processed_cores",
#                               core_cat[step]['x'][slct],core_cat[step]['y'][slct],core_cat[step]['z'][slct],
#                               core_cat[step]['vx'][slct],core_cat[step]['vy'][slct],core_cat[step]['vz'][slct],
#                               core_cat[step]['core_tag'][slct],
#                               core_cat[step]['infall_mass'][slct],
#                               core_cat[step]['infall_step'][slct])

                               
                               

    plt.figure()
    plt.hist(core_cat[step]['radius'][np.logical_not(np.isnan(core_cat[step]['radius']))],bins=np.logspace(-3,2,100))
    plt.xscale('log')
    plt.show()
    

    gal_r = core_cat[step]['radius']
    gal_m = core_cat[step]['infall_mass']
    slct = (gal_r < dis_len )& (gal_m > 10**infall_mass)
    gal_x = core_cat[step]['x'][slct]
    gal_y = core_cat[step]['y'][slct]
    gal_z = core_cat[step]['z'][slct]
    print("n2 start: ",gal_x.size)
    start = time.time()
    gal2_x, gal2_y, gal2_z, gal2_w, gal2_c= n2merger.n2merger3d(gal_x,gal_y,gal_z,merg_len,colors_out=True);
    end = time.time()
    print("halo merger time: ",end-start)
    start = time.time()
    #save_processed_core_cat(processed_core_loc,core_cat,slct,gal2_c)

    end = time.time()
    print("processed core catalog time: ",end-start)
    print("we got here")
    gal_xyz = return_xyz_formatted_array(gal2_x,gal2_y,gal2_z)
    r_bins = np.logspace(np.log10(0.001),np.log10(256.0/3.0),50)
    r_bins_avg = (r_bins[:1]+r_bins[:-1])/2.0
    r_val = tpcf(gal_xyz,r_bins,period=256.0,)
    print r_bins_avg.size,r_val.size
    start = time.time()
    #dtk_val,dtk_rad = dtk.autocorr3D(gal_x,gal_y,gal_z,128,256.0)
    end = time.time()
    print("dtk time: ",end-start)
    start = time.time()
    #dtk_n2_val,_,_ = dtk.autocorr3D_N2(gal_x[::1000],gal_y[::1000],gal_z[::1000],box_length=256.0,bins=r_bins)
    end = time.time()
    print("dtk time: ",end-start)
    #print dtk_val
    #print dtk_rad
    plt.figure()
    plt.plot(r_bins_avg,r_val,label='halo_tools')
    #plt.plot(dtk_rad,dtk_val,label='dtk')
    #plt.plot(r_bins_avg,dtk_n2_val,label='dtk n2')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('r [Mpc/h]')
    plt.ylabel('eta_gg')
    plt.grid()
    plt.legend()
    plt.show()

print("reading in gio files")
fof_cat.read_gio()
sod_cat.read_gio()

print("merging catalogs")
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
            print(i,"/",htags.size)
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


        max_r = cluster_radial_volume
        cores_slct = core_cat.cut_box(step,x,y,z,max_r,max_r,max_r,256,256,256,ignore_periodic=False)
        cores_x = core_cat[step]['x'][cores_slct]
        cores_y = core_cat[step]['y'][cores_slct]
        cores_z = core_cat[step]['z'][cores_slct]
        cores_r = core_cat[step]['radius'][cores_slct]
        cores_m = core_cat[step]['infall_mass'][cores_slct]

        #cutting sphere out:
        sphere_r2 = max_r**2
        cores_r2 = (cores_x-x)**2 + (cores_y-y)**2 + (cores_z-z)**2
        cores_slct = cores_r2<sphere_r2
        cores_x = cores_x[cores_slct]
        cores_y = cores_y[cores_slct]
        cores_z = cores_z[cores_slct]
        cores_m = cores_m[cores_slct]
        cores_r = cores_r[cores_slct]
        if(m200 > 1e15):
            indx = np.argmax(cores_m)
            print "making check", cores_m[indx],cores_r[indx]
        #cores2_x,cores2_y,cores2_z,cores2_r = core_cat.cut_box_kdtree(step,x,y,z,2)

        if(box):
            #indx = np.arange(cores_x.size)
            #numpy.random.shuffle(indx)
            #indx = np.argsort(cores_x)
            #cores_x = cores_x[indx].astype('float64')
            #cores_y = cores_y[indx].astype('float64')
            #cores_z = cores_z[indx].astype('float64')
            #cores_r = cores_r[indx].astype('float64')
            #cores_m = cores_m[indx].astype('float64')
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

        clstrs.append(clstr)


print("avg: ",np.average(avg))
print("avg: ",np.average(avg2))
print("done setting the damn things up")

zmr_core = zmr_valid.make_empty_zmr()

start_total = time.time()

start_time = time.time()

print("starting to calculate new zmr")
for i in range(0,len(clstrs)):

    clstr = clstrs[i]
    clstr.set_n2merg(n2merger)
    continue;
    dis_len = 0.3
    merg_len = 0.08
    infall_m = 11.61
    rad_cnt,rad_cnt_err,ngal = clstr.get_rad_bin_cnt(dis_len,merg_len)
    z_i, m_i = clstr.get_zi_mi()
    zmr_core['zmr_gal_cnt'][z_i,m_i,:] +=rad_cnt
    zmr_core['zm_counts'][z_i,m_i]+=1
    end_time = time.time()
    start_time = end_time

    #plotting individual 
    continue;
    res_x,res_y,res_z,res_w = clstr.process_cores_lib(dis_len,merg_len,infall_mass=infall_m,ndim_3d=True)
    print('before: ',clstr.core_x.size)
    print('after: ',res_x.size)
    rad_bin_area = zmr_valid.get_rad_bin_area()
    rad_bin_avg = zmr_valid.r_bin_avg
    rad_bin = zmr_valid.r_bins
    fig,ax = plt.subplots()
    #plt.title('m200=%.1e, zi=%d,mi=%d'%(clstr.m200,z_i,m_i))
    if(i%10==0):
        plt.show()
    plt.title('m200=%.1e'%(clstr.m200))
    plt.plot([],[],'k-',label='radial bins')

    plt.scatter(res_x,res_y,res_w*20.0,label='intact cores')
    plt.plot(clstr.x0,clstr.y0,'s')
    plt.plot(clstr.core_x,clstr.core_y,'g+',label='disrupted cores')
    plt.xlabel('x [mpc/h]')
    plt.ylabel('y [Mpc/h]')
    for r in rad_bin:
        c1 = plt.Circle((clstr.x0,clstr.y0),clstr.r200_comv*r,facecolor='None',edgecolor='k')
        ax.add_artist(c1)

    plt.axis('equal')
    plt.legend(loc='best',framealpha=0.5)

    continue
    plt.figure()
    plt.plot(rad_bin_avg,rad_cnt)
    plt.fill_between(rad_bin_avg,rad_cnt-rad_cnt_err,rad_cnt+rad_cnt_err,alpha='0.5')
    plt.gca().set_ylim(bottom=0)
    minimum = 0.1

    rad_den = rad_cnt/rad_bin_area+minimum
    rad_den_err = rad_cnt_err/rad_bin_area
    plt.figure()
    plt.plot(rad_bin_avg,rad_den)
    plt.fill_between(rad_bin_avg,np.clip(rad_den-rad_den_err,a_min=minimum,a_max=max(rad_den)),rad_den+rad_den_err,alpha='0.5')
    plt.yscale('log')
    plt.gca().set_ylim(bottom=minimum)

    if(i%10==0):
        plt.show()

plt.show()
def zmr_to_min(param):
    dis_len = param[0]
    merg_len = param[1]
    zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid)
    cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
    return cost

#taking infall mass as a variable
def zmr_to_min2(param):
    dis_len = param[0]
    merg_len = param[1]
    infall_m = param[2]
    zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid,infall_m)
    cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
    return cost



#set_min_var(clstrs,zmr_valid,zmr_sdss)

end_total = time.time()
print("total process: ",end_total-start_total)
start_total = time.time()
x0 = [.1,.05,11.60]
res = minimize(zmr_to_min2,x0,method='nelder-mead',
               options={'xtol':1e-5,'disp':True}) 

print("result: ",res)

dis_len = res.x[0]
merg_len = res.x[1]
infall_mass = res.x[2]
    
write_out_gal_clusters(clstrs,infall_mass,dis_len,merg_len);
#import corner
#print "Running MCMC now..."
#chain_len = 1000
#nwalker  = 100
#sampler = run_mcmc(res.x,nwalker,chain_len,zmr_sdss,zmr_valid,clstrs)

#samples = sampler.chain[:,chain_len/2:,:].reshape((-1,3))
#fig = corner.corner(samples, labels=['dis_len','merg_len','infall_m'])


end_total = time.time()
zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid,infall_mass)

print("fit time: ",end_total-start_total)


cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
print("this is the cost: ", cost)



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

#zmr_count_to_density(zmr_core,0)

colors = ['b','r','g','m','c','y']
c_i = 0
for i in range(0,zmr_core['z_bins'].size-1):
    if(np.sum(zmr_core['zm_counts'][i])==0):
        continue
    #galaxy density profile
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
c_i = 0;
#ngal density profile
zmr_mass_bins_avg = (zmr_core['mass_bins'][:-1]+zmr_core['mass_bins'][1:])/2.0
plt.figure()
for i in range(0,zmr_core['z_bins'].size-1):
    if(np.sum(zmr_core['zm_counts'][i])==0):
        continue;
    #cores
    mass_x = []
    ngal_y = []
    for j in range(0,zmr_core['mass_bins'].size-1):
        if(zmr_core['zm_counts'][i,j]>0):
            mass_x.append(zmr_mass_bins_avg[j])
            ngal_y.append(zmr_core['zm_ngal'][i,j])
    c = colors[c_i%len(colors)]
    c_i +=1
    plt.plot(mass_x,ngal_y,'-o',color=c,label='core %.2f<z<%.2f'%(zmr_core['z_bins'][i],zmr_core['z_bins'][i+1]))
    #cores
    mass_x = []
    ngal_y = []
    for j in range(0,zmr_core['mass_bins'].size-1):
        if(zmr_core['zm_counts'][i,j]>0):
            mass_x.append(zmr_mass_bins_avg[j])
            ngal_y.append(zmr_sdss['Ngal'][i,j])
    plt.plot(mass_x,ngal_y,'--x',color=c,label='sdss %.2f<z<%.2f'%(zmr_core['z_bins'][i],zmr_core['z_bins'][i+1]))
plt.xlabel('M200 [M_sun h^-1]')
plt.ylabel('Ngal')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.grid()
    
#make2dplot = raw_input("make 2d likelihood plot? [y/other]")

        




dtk.save_figs("figs/"+param.file+"/core_fit2.py/")

plt.show()
