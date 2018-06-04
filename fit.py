#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import dtk
import sod_fof_relation
from astropy.cosmology import WMAP9 as cosmo

param = dtk.Param(sys.argv[1])

zmr_sdss = np.load(param.get_string("zmr_file"))
for k in zmr_sdss.keys():
    print k

z_list   = zmr_sdss["z_bins"]
m_list   = zmr_sdss["mass_bins"]
r_list   = zmr_sdss["rad_bins"]
z_size   = z_list.size-1
m_size   = m_list.size-1
r_size   = r_list.size-1
z_avg    = (z_list[1:]+z_list[:-1])/2.0
m_avg    = (m_list[1:]+m_list[:-1])/2.0
r_avg    = (r_list[1:]+r_list[:-1])/2.0
core_file  = param.get_string("core_catalog_file")
core_steps = param.get_int_list("core_catalog_steps_load")
core_htag = []
core_x = []
core_y = []
core_z = []
core_size =[]
z_in = param.get_float("z_in")
z_out = param.get_float("z_out")
num_steps = param.get_int("num_steps")

a_in = 1.0/(z_in+1.0)
a_out = 1.0/(z_out+1.0)
a_del = (a_out-a_in)/(num_steps-1.0)
print a_in
print a_out
print a_del
fof_catalog_file = param.get_string("fof_catalog_file")
sod_catalog_file = param.get_string("sod_catalog_file")

def legal_z_index(zi):
    return zi >= 0 & zi<z_size

def legal_m_index(mi):
    return mi >= 0 & mi<m_size

def legal_r_index(ri):
    return ri >= 0 & ri<r_size

def load_core_data(core_file, step):
    print "Loading core step:", step
    fname = core_file.replace("${step}",str(step))
    core_htag.append(dtk.gio_read(fname,"fof_halo_tag"))
    core_x.append(dtk.gio_read(fname,"x"))
    core_y.append(dtk.gio_read(fname,"y"))
    core_z.append(dtk.gio_read(fname,"z"))
    core_size.append(dtk.gio_read(fname,"radius"))



class FOFCatalog:
    def __init__(self,fof_catalog_file,steps):
        self.steps = []
        self.fof_htags  = []
        self.fof_mass   = []
        for s in steps:
            self.load(fof_catalog_file,s)
    def load(self,fof_catalog_file,s):
        print "Loading fof catalog:", s
        fname = fof_catalog_file.replace("${step}",str(s))
        self.steps.append(s)
        htags = dtk.gio_read(fname,"fof_halo_tag")
        mass  = dtk.gio_read(fname,"fof_halo_mass")
        indx = np.argsort(htags)
        self.fof_htags.append(htags[indx])
        self.fof_mass.append(mass[indx])
    def lookup_mass(self,step,htag):
        step_indx = self.steps.index(step)
        indx = np.searchsorted(self.fof_htags[step_indx],htag)
        if(indx >= 0 and indx < self.fof_htags[step_indx].size):
            if(htag == self.fof_htags[step_indx][indx]):
                return self.fof_mass[step_indx][indx]
        return -1.0 # fof halo mass not found
    def mass_function(self):
        mass_bins = np.logspace(10,16,100)
        mass_bins_avg = (mass_bins[:-1]+mass_bins[1:])/2.0
        mass_counts  = np.zeros((len(self.steps),mass_bins_avg.size),dtype='i4')
        mass_counts2 = np.zeros((len(self.steps),mass_bins_avg.size),dtype='i4')
        hist         = np.zeros((len(self.steps),mass_bins_avg.size),dtype='i4')
        limit = 100000
        for i in range(0,len(self.steps)):
            print "working on %d/%d"%(i,len(self.steps))
            indx = np.digitize(self.fof_mass[i][:limit],mass_bins)-1
            indx2 = np.searchsorted(mass_bins,self.fof_mass[i][:limit])-1
            hist[i,:] = np.histogram(self.fof_mass[i][:limit],bins=mass_bins)[0]
            for j in range(0,self.fof_mass[i][:limit].size):
                if(indx[j]>=0 and indx[j]<mass_counts[i].size):
                    mass_counts [i,indx [j]]+=1
                if(indx2[j]>=0 and indx2[j]<mass_counts2[i].size):
                    mass_counts2[i,indx2[j]]+=1
        return mass_bins_avg,mass_counts,mass_counts2,hist

class SODCatalog:
    def __init__(self,sod_file,core_steps):
        self.core_steps = []
        self.fof_htags  = []
        self.sod_mass   = []
        self.sod_r200   = []
        for s in core_steps:
            self.load(sod_file,s)
    def load(self,sod_file,s):
        print "Loading sod catalog:", s
        fname = sod_file.replace("${num}",str(s))
        self.core_steps.append(s)
        htags = dtk.gio_read(fname,"fof_halo_tag")
        radius = dtk.gio_read(fname,"sod_radius")
        mass  = dtk.gio_read(fname,"sod_mass")
        indx = np.arg_sort(htags)
        self.htag.append(htags[indx])
        self.sod_mass.append(mass[indx])
        self.sod_radius.append(radius[indx])

    def lookup_radius(self,step,htag):
        if(not step in self.core_steps):
            return -1.0
        step_indx = self.steps.index(step)
        indx = np.searchsorted(fof_htags[step_indx],htag)
        if(indx >= 0 and indx < fof_htags[step_indx].size):
            if(htag == fof_htags[step_indx][indx]):
                return sod_radius[step_indx][indx]
        return -1.0 # fof halo mass not found
    def lookup_mass(self,step,htag):
        if(not step in self.core_steps):
            return -1.0
        step_indx = self.core_steps.index(step)
        indx = np.searchsorted(fof_htags[step_indx],htag) -1
        if(indx >= 0 and indx < fof_htags[step_indx].size):
            if(htag == fof_htags[step_indx][indx]):
                return sod_mass[step_indx][indx]
        return -1.0 # fof halo mass not found

        
def crit_density(z): #Msun/h /kpc^3
    gcm3_to_msunkpc3 = 1.477543e31
    density = cosmo.critical_density(z).value*gcm3_to_msunkpc3
    #print "crit desity(%f): %f Msun/kpc^3"%(z,density)
    #print "crit desity(%f): %e Msun/kpc^3"%(z,density*1000**3/cosmo.h**2)
    return density


def fof_mass_to_r200(fof_mass,z):
    return sod_mass_to_r200(fof_mass_to_sod_mass(fof_mass,z),z)

def sod_mass_to_r200(sod_mass,z):
    rho = crit_density(z)
    return ((sod_mass*3)/(4*200*np.pi*rho))**(1.0/3.0) 

def fof_mass_to_sod_mass(fof_mass,z):
    return sod_fof_relation.so_from_fof_z(fof_mass,z)

def group_by_id(data):
    # has to be sorted array
    if(data.size == 0):
        return np.zeros(0,dtype='i4'),np.zeros(0,dtype='i4')
    start = [0]
    val = data[0]
    for i in range(1,data.size):
        if(data[i] != val):
            start.append(i)
            val = data[i]
    group_start = np.array(start,dtype='i4')
    group_size = np.ones_like(group_start,dtype='i4')
    for i in range(0,group_start.size-1):
        group_size[i] = group_start[i+1]-group_start[i]
    group_size[-1]=data.size-group_start[-1]
    return group_start,group_size
        
def find_r200(step,fof_htag):
    rad = sod_cat.lookup_radius(fof_htag)
    if(rad ==-1.0):
        rad = fof_mass_to_r200(fof_cat.lookup_mass(step,fof_htag))
    return rad

def find_sod_mass(step,fof_htag):
    mass = sod_cat.lookup_mass(step,fof_htag)
    if(mass == -1.0):
        z = step_to_z(step)
        mass = fof_mass_to_sod_mass(fof_cat.lookup_mass(step,fof_htag),z)
    return mass

def step_to_z(step):
    a = a_in + step*a_del
    z = 1.0/a -1.0
    return z

def disrupt_merger_core_process(group_start,group_size,disruption_len,merger_len,step_indx):
    x = np.zeros(group_size)
    y = np.zeros(group_size)
    z = np.zeros(group_size)
    intact = np.zeros(group_size,dtype=bool)
    for i in range(0,group_size):
        x[i]=core_x[step_indx][group_start+i]
        y[i]=core_y[step_indx][group_start+i]
        z[i]=core_z[step_indx][group_start+i]
        intact[i]=core_size[step_indx][group_start+i]<disruption_len
    x = np.array(x[intact==1])
    y = np.array(y[intact==1])
    z = np.array(z[intact==1])
    if(x.size==0):
        return np.zeros(0),np.zeros(0),np.zeros(0)
    x,y,z,weight = merger_process_N2(x,y,z,merger_len)
    bcg_indx = np.argmax(weight)
    cen_x = x[bcg_indx]
    cen_y = y[bcg_indx]
    cen_z = z[bcg_indx]
    r_xy = np.sqrt((cen_x-x)**2+(cen_y-y)**2)
    r_xz = np.sqrt((cen_x-x)**2+(cen_z-z)**2)
    r_yz = np.sqrt((cen_y-y)**2+(cen_z-z)**2)
    return r_xy, r_xz, r_yz

def disrupt_merger_core_process2(group_start,group_size,disruption_len,merger_len,step_indx):
    #print group_start,group_size
    intact_slct = core_size[step_indx][group_start:group_start+group_size]<=disruption_len
    intact_size = np.sum(intact_slct)
    intact_indx = np.nonzero(intact_slct)
    #    print group_start,group_size,intact_slct,intact_size,intact_indx
    x = np.zeros(intact_size)
    y = np.zeros(intact_size)
    z = np.zeros(intact_size)
    j = 0
    for i in np.nonzero(intact_slct)[0]:
        x[j]=core_x[step_indx][group_start+i]
        y[j]=core_y[step_indx][group_start+i]
        z[j]=core_z[step_indx][group_start+i]
        j+=1
    if(x.size==0):
        return np.zeros(0),np.zeros(0),np.zeros(0)
    x,y,z,weight = merger_process_N2(x,y,z,merger_len)
    bcg_indx = np.argmax(weight)
    cen_x = x[bcg_indx]
    cen_y = y[bcg_indx]
    cen_z = z[bcg_indx]
    r_xy = np.sqrt((cen_x-x)**2+(cen_y-y)**2)
    r_xz = np.sqrt((cen_x-x)**2+(cen_z-z)**2)
    r_yz = np.sqrt((cen_y-y)**2+(cen_z-z)**2)
    return r_xy, r_xz, r_yz

def merger_process_N2(x,y,z,merger_len):
    if(x.size==0):
        return np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0,dtype='i4')
    if(merger_len==0):
        return x,y,z,np.ones_like(x)
    to_plot = x.size>20000
    if(to_plot):
        plt.figure()
        plt.scatter(x,y,s=80,c='b',marker='x',label='input')
    color = np.arange(0,x.size)
    merg_len = merger_len**2
    for i in range(0,x.size):
        for j in range(0,i):
            dist = (x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2
            if(dist <= merg_len and color[i] != color[j]):
                merge_colors(color,color[i],color[j])
    indx = np.argsort(color)
    color = color[indx]
    x = x[indx]
    y = y[indx]
    z = z[indx]
    move_together(x,y,z)
    [group_start,group_size] = group_by_id(color)
    res_x = np.zeros(group_start.size)
    res_y = np.zeros(group_start.size)
    res_z = np.zeros(group_start.size)
    res_weight = np.zeros(group_start.size)
    for i in range(0,group_start.size):
        res_x[i] = np.average(x[group_start[i]:group_start[i]+group_size[i]])
        res_y[i] = np.average(y[group_start[i]:group_start[i]+group_size[i]])
        res_z[i] = np.average(z[group_start[i]:group_start[i]+group_size[i]])
        res_weight[i] = group_size[i]
    if(to_plot):
        plt.scatter(res_x,res_y,s=80,c='None',edgecolors='r',label='result')
        plt.title("merger_len: %.3f\nbefore:%d after:%d"%(merger_len,x.size,res_x.size))
        plt.show()
    return res_x,res_y,res_z,res_weight
            
def merge_colors(color,c1,c2):
    color[color==c2]=c1

def move_together(x,y,z):
    if(x.size == 0):
        return
    x1 = x[0]
    y1 = y[0]
    z1 = z[0]
    x_diff = x-x1
    y_diff = y-y1
    z_diff = z-z1
    for i in range(0,x.size):
        if(x_diff[i] > 128.0):
            x[i] -= 256.0
        elif(x_diff[i] < -128.0):
            x[i] += 256.0
        if(y_diff[i] > 128.0):
            y[i] -= 256.0
        elif(y_diff[i] < -128.0):
            y[i] += 256.0
        if(z_diff[i] > 128.0):
            z[i] -= 256.0
        elif(z_diff[i] < -128.0):
            z[i] += 256.0

def populate_mr(mr_gal_cnt,m_cnts,rad,sod_mass):
    m_indx = np.searchsorted(m_list,sod_mass)-1
    Ngal = 0
    r_gal = np.zeros((r_size),dtype='i4')
    if(m_indx >= 0 and m_indx < m_size):
        m_cnts[m_indx] += 1.0
        rad_indxes = np.searchsorted(r_list,rad)-1
        for r in rad_indxes:
            if(r>=0 and r<r_size):
                mr_gal_cnt[m_indx,r] +=1
                Ngal +=1
                r_gal[r] +=1
    return Ngal,r_gal,m_indx

def print_z_indx():
    for step_indx in range(0,len(core_steps)):
        z = step_to_z(core_steps[step_indx])
        z_indx = np.searchsorted(z_list,z)-1
        print "z:",z,"->",z_indx

def make_core_zmr(disruption_len, merger_len):
    print "Making zmr",z_size,m_size,r_size
    zmr_make_timer = time.time()
    zm_Ngal       = np.zeros((z_size,m_size))
    zm_Ngal_err   = np.zeros((z_size,m_size))
    zm_Ngal_var   = np.zeros((z_size,m_size))
    zmr_gal_count = np.zeros((z_size,m_size,r_size),dtype='i4')
    zmr_gal_density = np.zeros((z_size,m_size,r_size),dtype='f4')
    zmr_gal_density_var = np.zeros((z_size,m_size,r_size),dtype='f4')
    zmr_gal_density_err = np.zeros((z_size,m_size,r_size),dtype='f4')
    zm_Ngal_var_list = {}
    zmr_gal_var_list = {}
    for i in range(0,z_size):
        for j in range(0,m_size):
            zm_Ngal_var_list[i,j]=[];
            for k in range(0,r_size):
                zmr_gal_var_list[i,j,k] = [];
    zm_counts = np.zeros((z_size,m_size),dtype='i4')
    for step_indx in range(0,len(core_steps)):
        print "starting..."
        z = step_to_z(core_steps[step_indx])
        print "wtf"
        z_indx = np.searchsorted(z_list,z)-1
        print z_indx,z_list,z
        print "Getting host halos..."
        group_start,group_size = group_by_id(core_htag[step_indx])
        print "done. avarage size: ",np.average(group_size)
        print "average n2: ", np.average(group_size**2)
        mr_gal_counts  = np.zeros((m_size,r_size),dtype='i4')
        m_counts      = np.zeros((m_size),dtype='i4')
        start = time.time()
        for j in range(0,group_start.size):
            if(j%10000==9999):
                end = time.time()
                print "%d/%d"%(j+1,group_start.size), (end-start)
                start = end

            host_htag = core_htag[step_indx][group_start[j]]
            sod_mass = find_sod_mass(core_steps[step_indx],host_htag)
            if(sod_mass<5e13):
                continue;
            rad_xy,rad_xz,rad_yz = disrupt_merger_core_process2(group_start[j],group_size[j],disruption_len,merger_len,step_indx)
            Ngal,r_gal,m_indx = populate_mr(mr_gal_counts,m_counts,rad_xy,sod_mass)
            if(legal_m_index(m_indx)):
                zm_Ngal_var_list[z_indx,m_indx].append(Ngal)
                for k in range(0,r_size):
                    zmr_gal_var_list[z_indx,m_indx,k].append(r_gal[k])
        zmr_gal_count[z_indx,:,:] += mr_gal_counts
        zm_counts[z_indx,:] += m_counts
    zmr_gal_count = zmr_gal_count.astype('f4')
    zm_Ngal = np.nansum(zmr_gal_count,axis=2)/zm_counts
    radial_bin_area = np.pi*(r_list[1:]**2-r_list[:-1]**2)
    for i in range(0,r_size):
        zmr_gal_count [:,:,i] = zmr_gal_count[:,:,i]/zm_counts
        zmr_gal_density[:,:,i] = zmr_gal_count[:,:,i]/radial_bin_area[i]
    for i in range(0,z_size):
        for j in range(0,m_size):
            zm_Ngal_var[i,j]=np.std(zm_Ngal_var_list[i,j])
            zm_Ngal_err[i,j]=zm_Ngal_var[i,j]/np.sqrt(len(zm_Ngal_var_list[i,j]))
            for k in range(0,r_size):
                zmr_gal_density_var[i,j,k]=np.std(zmr_gal_var_list[i,j,k])/radial_bin_area[k]
                zmr_gal_density_err[i,j,k]=zmr_gal_density_var[i,j,k]/np.sqrt(zmr_gal_density_var[i,j,k])
    zmr = {}
    zmr['zm_Ngal'] = zm_Ngal
    zmr['zm_Ngal_var'] = zm_Ngal_var
    zmr['zm_Ngal_err'] = zm_Ngal_err
    zmr['zmr_gal_density'] = zmr_gal_density
    zmr['zmr_gal_density_var'] = zmr_gal_density_var
    zmr['zmr_gal_density_err'] = zmr_gal_density_err
    zmr['zm_counts']  =zm_counts
    print "time for single val: ", time.time() - zmr_make_timer
    return zmr

def calc_Ngal_cost(core_zmr,sdss_zmr):
    #just doing Ngal for the moment
    result = 0.0
    for i in range(0,z_size):
        for j in range(0,m_size):
            err = np.sqrt(core_zmr['zm_Ngal_err'][i,j]**2 + sdss_zmr['Ngal_err'][i,j]**2)
            diff = np.abs(core_zmr['zm_Ngal'][i,j]-sdss_zmr['Ngal'][i,j])
            cost = (diff/err)**2
            if(core_zmr['zm_counts'][i,j]==0 or sdss_zmr['zmr_cnt'][i,j,0]==0):
                cost = 0.0
            result += cost
    return result

def calc_gal_density_cost(core_zmr,sdss_zmr):
    result =0.0
    for i in range(0,z_size):
        for j in range(0,m_size):
            for k in range(0,r_size):
                diff = np.abs(core_zmr['zmr_gal_density'][i,j,k] - sdss_zmr['rad_prof'][i,j,k])
                err  = core_zmr['zmr_gal_density'][i,j,k]**2 + sdss_zmr['rad_prof_err'][i,j,k]**2
                cost = diff**2/err
                if(core_zmr['zm_counts'][i,j] == 0 or sdss_zmr['zmr_cnt'][i,j,k] == 0):
                    cost = 0.0
                result += cost
    return result

print "Printing z index"
print_z_indx()

print core_steps
fof_cat = FOFCatalog(fof_catalog_file,core_steps)
sod_cat = SODCatalog(sod_catalog_file,[])



for s in core_steps:
    load_core_data(core_file,s)
    
# for i in range(0,len(core_steps)):
#     htags = np.unique(core_htag[i])
#     fof_masses = []
#     sod_masses = []
#     for h in htags:
#         fof_masses.append(fof_cat.lookup_mass(core_steps[i],h))
#         sod_masses.append(find_sod_mass(core_steps[i],h))
#     plt.figure()
#     plt.plot(m_avg,np.histogram(fof_masses,bins=m_list)[0])
#     plt.plot(m_avg,np.histogram(sod_masses,bins=m_list)[0])
#     plt.yscale('log')
#     plt.xscale('log')
#     not_found = np.sum(fof_masses==-1.0)
#     found     = np.sum(fof_masses!=-1.0)
#     plt.title("Found: %d, not found: %d"%(found,not_found))
# plt.show()
dis_bins = np.logspace(np.log10(0.03),np.log10(0.4),5)
dis_avg  = np.sqrt(dis_bins[1:]*dis_bins[:-1])
merg_bins = np.logspace(np.log10(0.001),np.log10(0.3),5)
merg_avg = np.sqrt(merg_bins[1:]*merg_bins[:-1])


def zmr_diff(param):
    disruption_len = param[0]
    merger_len  = param[1]
    print "\n\ntrying dis: ",disruption_len," merger_len: ",merger_len
    core_zmr = make_core_zmr(disruption_len,merger_len)
    cost_mat_gal_den = calc_gal_density_cost(core_zmr,zmr_sdss)
    print "cost: ", cost_mat_gal_den
    return cost_mat_gal_den

from scipy.optimize import minimize

def plot_zmr(zmr,zmr_err,zmr_var,z_is,m_is):
    for z_i in z_is:
        for m_i in m_is:
            plt.plot(r_bins_avg,zmr['zmr_gal_density'][z_i,m_i],label='z=%.1f,m=%.1e'%(z_bins_avg[z_i],m_bins_avg[m_i]))
            plt.yerror(r_bins_avg,zmr['zmr_gal_density_err'][z_i,m_i],(z_bins_avg[z_i],m_bins_avg[m_i]))
    return 


dis_len  = .2
merg_len = .2
colors = ['b','r','g','m','y']
zmr_fit = make_core_zmr(dis_len,merg_len)
print zmr_fit.keys()
print zmr_sdss.keys()
print "\n"
for z_i in range(0,len(z_avg)):
    has_fig = False
    print z_i
    for m_i in range(0,len(m_avg)):
        color = colors[m_i%len(colors)]
        if(np.sum(zmr_fit['zmr_gal_density'][z_i,m_i]>0)>0):
            if(np.sum(zmr_sdss['rad_prof'][z_i,m_i]>0)>0):
                if(not has_fig):
                        plt.figure()
                        plt.title('z=%.1f'%z_avg[z_i])
                        has_fig = True
                plt.plot(r_avg,zmr_fit['zmr_gal_density'][z_i,m_i],'o-',label='m=%.1e'%m_avg[m_i],c=color)
                plt.errorbar(r_avg,zmr_fit['zmr_gal_density'][z_i,m_i],yerr=zmr_fit['zmr_gal_density_err'][z_i,m_i],c=color)
                plt.plot(r_avg,zmr_sdss['rad_prof'][z_i,m_i],'x--',label='m=%.1e'%m_avg[m_i],c=color)
                plt.errorbar(r_avg,zmr_sdss['rad_prof'][z_i,m_i],yerr=zmr_sdss['rad_prof_err'][z_i,m_i],c=color)
    if(has_fig):
        plt.yscale('log')
        plt.ylim([0,200])
        plt.legend(loc='best')
plt.show()
exit()
print "starting to minimize the params!"
x0 = np.array([.2,.2])
res = minimize(zmr_diff,x0,method='nelder-mead',
               options={'xtol':1e-2,'disp':True})
print res
zmr_fit = make_core_zmr(res[0],res[1])

for z_i in range(1,2):
    plt.figure()
    plt.title('z=%.1f'%z_bins_avg[z_i])
    print z_i
    for m_i in range(0,len(m_bins_avg)):
        if(np.sum(zmr_fit['zmr_gal_density'][z_i,m_i]>0)>0):
            plt.plot(r_bins_avg,zmr_fit['zmr_gal_density'][z_i,m_i],'o',label='m=%.1e')
        if(np.sum(zmr_sdss['zmr_gal_desity'][z_i,m_i]>0)>0):
            plt.plot(r_bins_avg,zmr_sdss['zmr_gal_desity'][z_i,m_i],'x',label='m=%.1e')
    plt.yscale('log')
    plt.ylim([0,100])

exit()
cost_mat_Ngal = np.zeros((dis_avg.size,merg_avg.size))
cost_mat_gal_den = np.zeros((dis_avg.size,merg_avg.size))
print "starting the thing"
for i in range(0,dis_avg.size):
    for j in range(0,merg_avg.size):
        print "\n\ni: %d/%d, j: %d/%d"%(i+1,dis_avg.size,j+1,merg_avg.size)
        core_zmr = make_core_zmr(dis_avg[i],merg_avg[j])
        cost_mat_Ngal[i,j] = calc_Ngal_cost(core_zmr,zmr_sdss)
        cost_mat_gal_den[i,j] = calc_gal_density_cost(core_zmr,zmr_sdss)
print calc_Ngal_cost(core_zmr,zmr_sdss)
print "Done."
from matplotlib.colors import LogNorm

plt.figure()
X,Y = np.meshgrid(dis_bins,merg_bins)
plt.title("Ngal Cost")
plt.pcolor(X,Y,cost_mat_Ngal.T,vmin=np.min(cost_mat_Ngal),vmax=np.max(cost_mat_Ngal),norm=LogNorm())
plt.xlabel('disruption length')
plt.ylabel('merger length')
plt.yscale('log')
plt.xscale('log')
plt.colorbar()

plt.figure()
plt.title("Rad prof Cost")
plt.pcolor(X,Y,cost_mat_gal_den.T,norm=LogNorm(),vmax = min(1e2,np.max(cost_mat_gal_den)))
plt.xlabel('disruption length')
plt.ylabel('merger length')
plt.yscale('log')
plt.xscale('log')
plt.colorbar()

plt.figure()
plt.title("Multi of cost")
plt.pcolor(X,Y,np.sqrt(cost_mat_gal_den.T*cost_mat_Ngal.T),norm=LogNorm(),vmax = min(1e2,np.max(cost_mat_gal_den)))
plt.xlabel('disruption length')
plt.ylabel('merger length')
plt.yscale('log')
plt.xscale('log')
#plt.colorbar()

plt.figure()
plt.title("Sum of cost")
plt.pcolor(X,Y,(cost_mat_gal_den.T+cost_mat_Ngal.T)/2.0,norm=LogNorm(),vmax = min(1e2,np.max(cost_mat_gal_den)))
plt.xlabel('disruption length')
plt.ylabel('merger length')
plt.yscale('log')
plt.xscale('log')
#plt.colorbar()

dtk.save_figs("figs/"+param.file+"/")
plt.show()

exit()
plt.figure()
for i in range(0,7):#z_size):
    plt.plot(m_avg,core_zmr['zm_Ngal'][i,:],'-x',label='z:%.2f'%(z_avg[i]))
    plt.plot(m_avg,zmr_sdss['Ngal'][i,:],"-s",label='z:%.2f'%(z_avg[i]))
plt.title("Ngal")
plt.legend()    
plt.yscale('log')
plt.xscale('log')

plt.show()
plt.figure()
for i in range(0,7):
    plt.plot(m_avg,core_zmr['zm_counts'][i,:],'-x',label='z:%.2f'%(z_avg[i]))
    plt.plot(m_avg,zmr_sdss['zm_cnt'][i,:],'-s',label='z:%.2f'%(z_avg[i]))
plt.title("halo mass count")
plt.legend()
plt.yscale('log')
plt.xscale('log')


plt.show()
