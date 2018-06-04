#!/usr/bin/env python2.7
import numpy as np
import scipy as sp
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from core_fit2_util import *
import dtk

def frag_to_real(htags):
    return np.abs(htags) & 0x0000ffffffffffff


class Catalog:
    def __init__(self,file_loc=None,step_string=None):
        self.step_data = {}
        self.data_names = []
        self.file_source_ = ""
        self.step_string_ = ""
        self.srt = None
        if(file_loc != None):
            if(step_string != None):
                self.set_file(file_loc,step_string)
            else:
                self.set_file(file_loc)


    def set_file(self,file_source,step_string="${step}"):
        self.file_source_ = file_source
        self.step_string_ = step_string

    def add_steps(self,steps):
        for step in steps:
            self.step_data[step]= {}

    def add_step(self,step):
        self.step_data[step]={}

    def get_steps(self):
        return self.step_data.keys()
        
    def add_data_name(self, name):
        self.data_names.append(name)

    def read_gio(self,verbose=False):
        if(verbose):
            print "from file",self.file_source_
        for step in self.step_data.keys():
            if(verbose):
                print "\treading in step: ",step
            file_name = self.file_source_.replace(self.step_string_,"%d"%step)
            for name in self.data_names:
                if(verbose):
                    print "\treading in ",name
                data = dtk.gio_read(file_name,name)
                self.step_data[step][name] = data

    def read_hdf5(self,verbose=False):
        print "from file",self.file_source_
        for step in self.step_data.keys():
            print "\treading in step: ",step
            file_name = self.file_source_.replace(self.step_string_,"%d"%step)
            for name in self.data_names:
                print "\treading in ",name
                print "not implmented!!!!"
    
    def __getitem__(self,key):
        return self.step_data[key]

    def join(self,catalog1,catalog2,join_on='fof_halo_tag',verbose=False):
        #only merger if they are on the same timestep
        if(verbose):
            print "merging catalogs"

        for step in catalog1.step_data.keys():
            if step in catalog2.step_data.keys():
                if(verbose):
                    print "Both have step: ",step
                self.step_data[step]= {}
                srt1 = np.argsort(catalog1.step_data[step][join_on])
                srt2 = np.argsort(catalog2.step_data[step][join_on])
                if(verbose):
                    print "cat1 size: ", srt1.size
                    print "cat2 size: ", srt2.size
                i1 = 0
                i2 = 0
                i1_max = srt1.size
                i2_max = srt2.size
                match1 =[]
                match2 =[]
                unmatched1 = []
                unmatched2 = []
                if(verbose):
                    print "sorting"
                while(i1<i1_max and i2<i2_max):
                    val1 = catalog1.step_data[step][join_on][srt1[i1]]
                    val2 = catalog2.step_data[step][join_on][srt2[i2]]
                    if(val1 == val2):
                        match1.append(srt1[i1])
                        match2.append(srt2[i2])
                        i1+=1
                        i2+=1
                    elif(val1 > val2):
                        unmatched2.append(srt2[i2])
                        i2+=1
                    else: #(val1 < val2):
                        unmatched1.append(srt1[i1])
                        i1+=1
                while(i1<i1_max):
                    unmatched1.append(srt1[i1])
                    i1+=1
                while(i2<i2_max):
                    unmatched2.append(srt2[i2])
                    i2+=1
                match1 = np.array(match1)
                match2 = np.array(match2)
                unmatched1 = np.array(unmatched1)
                unmatched2 = np.array(unmatched2)
                if(verbose):
                    print "match1: ", match1.size
                    print "match2: ", match2.size
                    print "unmatched1: ", unmatched1.size
                    print "unmatched2: ", unmatched2.size

                    print "making merged catalog"

                #copy over the matched fields to this catalog
                #and only leave the unmatched data rows in the
                #original catalog
                names_this_step = []
                for name in catalog1.data_names:
                    self.step_data[step][name] = catalog1[step][name][match1]
                    if(verbose):
                        print step,name,"size:",self.step_data[step][name].size
                    names_this_step.append(name)
                    if(name not in self.data_names):
                        self.data_names.append(name)
                    if(unmatched1.size != 0):
                        catalog1[step][name]=catalog1[step][name][unmatched1]
                    else:
                        catalog1[step][name]=np.zeros(0,dtype=catalog1[step][name].dtype)
                for name in catalog2.data_names:
                    print "name:",name
                    if(name not in names_this_step): #not to add the same column multiple times
                        self.step_data[step][name] = catalog2[step][name][match2]
                    if(name not in self.data_names):
                        self.data_names.append(name)
                    if(unmatched2.size != 0):
                        catalog2[step][name]= catalog2[step][name][unmatched2]
                    else:
                        catalog2[step][name]=np.zeros(0,dtype=catalog2[step][name].dtype)
                if(verbose):
                    print "step: ", step, "vars: ", self.step_data[step].keys()
        if(verbose):
            print "\n\n"
            for step in self.step_data.keys():
                print "step: ",step,"vars: ", self.step_data[step].keys()

    def apply_function(self,var_name,function,*args):
        for step in self.step_data.keys():
            self.step_data[step][var_name] = function(self.step_data[step][var_name],*args)

    def sort(self,sort_by):
        self.srt = {}
        self.srt_var = sort_by
        for step in self.step_data.keys():
            self.srt[step] = np.argsort(self.step_data[step][sort_by])
        
    def find(self,step,val):
        if(self.srt==None):
            print "Not sorted yet"
            raise Exception('Not Sorted yet!')
        srt_indx = np.searchsorted(self.step_data[step][self.srt_var],val,sorter=self.srt[step])
        if(srt_indx >=0 and srt_indx < self.srt[step].size):
            indx = self.srt[step][srt_indx]
            if(self.step_data[step][self.srt_var][indx] == val):
                return indx
            else:
                return np.nan
        else:
            return np.nan
    
    def find_all(self,step,val):
        if(self.srt==None):
            print "Not sorted yet"
            raise Exception('Not Sorted yet!')
        srt_indx_start = np.searchsorted(self.step_data[step][self.srt_var],val,sorter=self.srt[step],side='left')
        srt_indx_end = np.searchsorted(self.step_data[step][self.srt_var],val,sorter=self.srt[step],side='right')-1
        if(np.isfinite(srt_indx_start) and np.isfinite(srt_indx_end)):
            val1 = self.step_data[step][self.srt_var][self.srt[step][srt_indx_start]]
            val2 = self.step_data[step][self.srt_var][self.srt[step][srt_indx_end]]
            if(val1 == val and val2 == val):
                result = []
                for srt_indx in range(srt_indx_start,srt_indx_end+1,1):
                    result.append(self.srt[step][srt_indx])
                return np.array(result)
            else:
                return np.atleast_1d(np.array((),dtype=int))
        return  np.atleast_1d(np.array((),dtype=int))

    def cut_box(self,step,x0,y0,z0,x_lim,y_lim,z_lim,x_wrap,y_wrap,z_wrap,var_x='x',var_y='y',var_z='z',ignore_periodic=False):
        slct_x = self.cut_dim(step,x0,x_lim,x_wrap,var_x=var_x,ignore_periodic=ignore_periodic)
        slct_y = self.cut_dim(step,y0,y_lim,y_wrap,var_x=var_y,ignore_periodic=ignore_periodic)
        slct_z = self.cut_dim(step,z0,z_lim,z_wrap,var_x=var_z,ignore_periodic=ignore_periodic)
        slct= slct_x & slct_y & slct_z
        return slct

    def cut_dim(self,step,x0,x_lim,x_wrap,var_x='x',ignore_periodic=False):
        if(ignore_periodic):
            x_shifted = self.step_data[step][var_x]
        else:
            x_shifted=self.shift_x(x0,self.step_data[step][var_x],x_wrap)
        dist = x_shifted - x0
        slct = (dist<x_lim) & (dist>-x_lim)
        return slct

    def shift_x(self,x0,x,x_wrap):
        x1 = np.copy(x)
        x_pos = x1+x_wrap
        x_neg = x1-x_wrap
        dist = x-x0
        result = np.where(dist>x_wrap/2.0,x_neg,x1)
        result2 = np.where(dist<-x_wrap/2.0,x_pos,result)
        return result2
        
    def make_kdtree(self,var_x='x',var_y='y',var_z='z'):
        self.kd_data = {}
        self.kdtree = {}
        for step in self.step_data.keys():
            self.kd_data[step]  = np.zeros((self.step_data[step][var_x].size,3))
            self.kd_data[step][:,0] = self.step_data[step][var_x]
            self.kd_data[step][:,1] = self.step_data[step][var_y]
            self.kd_data[step][:,2] = self.step_data[step][var_z]
            self.kdtree[step]       = sp.spatial.KDTree(self.kd_data[step])

    def cut_box_kdtree(self,step,x0,y0,z0,r0):
        #the box diagonal is a factor sqrt(3) longer than the length
        indxs = self.kdtree[step].query_ball_point((x0,y0,z0),r0*np.sqrt(3))
        #print "cut box indxs: ",indxs
        x = self.kd_data[step][indxs,0]
        y = self.kd_data[step][indxs,1]
        z = self.kd_data[step][indxs,2]
        r = self.step_data[step]['radius'][indxs]
        slct_x = ((x-x0)<r0) & ((x-x0)>-r0)
        slct_y = ((y-y0)<r0) & ((y-y0)>-r0)
        slct_z = ((z-z0)<r0) & ((z-z0)>-r0)
        slct = slct_x & slct_y & slct_z
        return x[slct],y[slct],z[slct],r[slct]

def frag_to_real(htags):
    return np.abs(htags) & 0x0000ffffffffffff


if __name__=="__main__":
    fof_loc = "/media/luna1/dkorytov/data/AlphaQ/fof/m000-${step}.fofproperties"
    sod_loc = "/media/luna1/dkorytov/data/AlphaQ/sod/m000-${step}.sodproperties"
    core_loc = "/media/luna1/dkorytov/data/AlphaQ/core_catalog/12_20_16.AlphaQ.${step}.coreproperties"
    fof_cat = Catalog()
    sod_cat = Catalog()
    merg_cat = Catalog()
    core_cat = Catalog(core_loc)
    fof_cat.set_file(fof_loc,"${step}")
    sod_cat.set_file(sod_loc,"${step}")
    #fof_cat.add_steps((499))#,453,401))
    #sod_cat.add_steps((499))#,453,401))
    step1 = 102
    fof_cat.add_step(step1)
    sod_cat.add_step(step1)
    core_cat.add_step(step1)
    fof_cat.add_data_name("fof_halo_tag")
    fof_cat.add_data_name("fof_halo_mass")
    #fof_cat.add_data_name("fof_halo_center_x")
    #fof_cat.add_data_name("fof_halo_center_y")
    #fof_cat.add_data_name("fof_halo_center_z")

    sod_cat.add_data_name("fof_halo_tag")
#    sod_cat.add_data_name("fof_halo_count")
    sod_cat.add_data_name("sod_halo_mass")
    #sod_cat.add_data_name("sod_halo_min_pot_x")
    #sod_cat.add_data_name("sod_halo_min_pot_y")
    #sod_cat.add_data_name("sod_halo_min_pot_z")
    core_cat.add_data_name("fof_halo_tag")
    core_cat.add_data_name("x")
    core_cat.add_data_name("y")
    core_cat.add_data_name("z")
    core_cat.add_data_name("radius")

    fof_cat.read_gio()
    sod_cat.read_gio()
    print "\n\n reading in core catalog"
    core_cat.read_gio()
    print "sorting core_catalog"
        
    val = 747856
    core_cat.apply_function('fof_halo_tag',frag_to_real)
    core_cat.sort("fof_halo_tag")
    core_indx = core_cat.find_all(step1,val)
    core_slct = core_cat[step1]['fof_halo_tag']==val
    print np.sum(core_slct)
    print core_cat[step1]['fof_halo_tag'][core_indx]
    print core_cat[step1]['x'][core_indx]
    print core_cat[step1]['y'][core_indx]
    print core_cat[step1]['z'][core_indx]
    cx = 128.0
    cy = 128.0
    cz = 128.0
    slct_vol = core_cat.cut_box(step1,cx,cy,cz,20.0,20.0,20.0,256.0,256.0,256.0)
    xs = core_cat[step1]['x'][slct_vol]
    ys = core_cat[step1]['y'][slct_vol]
    zs = core_cat[step1]['z'][slct_vol]
    print np.sum(slct_vol)
    plt.figure()
    plt.plot(xs,ys,'x')
    plt.plot(cx,cy,'s')
    slct_z = (core_cat[step1]['z']<20.0)|(core_cat[step1]['z']>(256.0-20.0))
    slct_z = (core_cat[step1]['z']>(cz-20.0))&(core_cat[step1]['z']<(cz+20))
    plt.plot(core_cat[step1]['x'][slct_z],core_cat[step1]['y'][slct_z],'+')
    #kdtrees
    core_cat.make_kdtree()
    x2,y2,z2,r2 = core_cat.cut_box_kdtree(step1,cx,cy,cz,20)
    plt.plot(x2,y2,'x')
    plt.grid()
    plt.show()
    exit()
    zmr = {}
    zmr['z_bins'] = np.linspace(0,2,10)
    zmr['mass_bins'] = np.logspace(1e11,2e15,10)
    zmr['rad_bins']  = np.linspace(0,10,10)
    print "Making zmr"
    zmr_validator = ZMRIndexValidator(zmr)
    print "making cluster objct"
    clstr = Cluster(0.1,1e13,1,cx,cy,cz,zmr_validator)
    print "setting cores"
    clstr.set_cores(core_cat[step1]['x'][slct_vol],
                    core_cat[step1]['y'][slct_vol],
                    core_cat[step1]['z'][slct_vol],
                    core_cat[step1]['radius'][slct_vol])
    print "done setting cores"
    res_x, res_y, res_weight = clstr.process_cores(1.0,1.101)
    print clstr.core_x
    print clstr.core_y
    print res_x
    print res_y
    plt.scatter(res_x,res_y,s=res_weight*20,c='c')




    merg_cat.join(fof_cat,sod_cat,verbose=False)
    print merg_cat[step1].keys()
    print "print catalog data var"
    for names in merg_cat.data_names:
        print names,
    print ""
    for i in range(0,5):
        for name in merg_cat.data_names:
            print merg_cat[step1][name][i],
        print ""
    
    print "sorting"
    merg_cat.sort("fof_halo_tag")
    val = 747856
    indx = merg_cat.find(step1,val)
    print "finding val ", indx
    print indx == np.nan
    print np.isnan(indx)
    print merg_cat[step1]['fof_halo_tag'][indx]
    mass_bins = np.logspace(11,15,100)
    x_0 = 0.226259
    y_0 = 186.099 
    z_0 = 86.4152
    exit()
    plt.figure()    
    merg_cat.apply_function("sod_halo_mass",lambda x:(x*100.0))
    print "this is the result:"
    print merg_cat[step1]["fof_halo_mass"][:10]
    print merg_cat[step1]["sod_halo_mass"][:10]
    print merg_cat[step1]["sod_halo_mass"][:10]/merg_cat[step1]["fof_halo_mass"][:10]
    H,_,_ = np.histogram2d(merg_cat[step1]["fof_halo_mass"],merg_cat[step1]["sod_halo_mass"],bins=(mass_bins,mass_bins))
    plt.pcolormesh(mass_bins,mass_bins,H.T,cmap='PuBu',norm=colors.LogNorm())
    plt.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],'r')
    plt.xlabel('FoF Halo Mass [Msun/h]')
    plt.ylabel('SOD Halo Mass [Msun/h]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()

    ratio_bins = np.linspace(0,2,100)
    plt.figure()
    ratio = merg_cat[step1]["fof_halo_mass"]/merg_cat[step1]["sod_halo_mass"]
    H,_,_ = np.histogram2d(merg_cat[step1]["fof_halo_mass"],ratio,bins=(mass_bins,ratio_bins))
    plt.pcolormesh(mass_bins,ratio_bins,H.T,cmap='PuBu',norm=colors.LogNorm())
    plt.xlabel('FoF Halo Mass [Msun/h]')
    plt.ylabel('fof/sod mass')
    plt.xscale('log')
    plt.grid()
    
    plt.show()
                           
