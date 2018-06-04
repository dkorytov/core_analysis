
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import dtk
import scipy.spatial
from scipy.optimize import minimize
import ctypes as ct
import os
import h5py

#Define where the library is and load it
class N2Merger:
    def __init__(self,loc):
        self.lib = ct.CDLL(os.path.abspath(loc))
        #setup the return and argument types
        self.lib.n2_merger_double.restype=None
        self.lib.n2_merger_double.argtypes=[ct.POINTER(ct.c_double),
                                            ct.POINTER(ct.c_double),
                                            ct.POINTER(ct.c_int),
                                            ct.POINTER(ct.c_int),
                                            ct.c_double,
                                            ct.POINTER(ct.c_int64)]
        self.lib.n2_merger_float.restype=None
        self.lib.n2_merger_float.argtypes=[ct.POINTER(ct.c_float),
                                           ct.POINTER(ct.c_float),
                                           ct.POINTER(ct.c_int),
                                           ct.POINTER(ct.c_int),
                                           ct.c_float,
                                           ct.POINTER(ct.c_int64)]
        self.lib.n2_merger_double3d.restype=None
        self.lib.n2_merger_double3d.argtypes=[ct.POINTER(ct.c_double),
                                              ct.POINTER(ct.c_double),
                                              ct.POINTER(ct.c_double),
                                              ct.POINTER(ct.c_int32),
                                              ct.POINTER(ct.c_int32),
                                              ct.c_double,
                                              ct.POINTER(ct.c_int64)]
        self.lib.n2_merger_float3d.restype=None
        self.lib.n2_merger_float3d.argtypes=[ct.POINTER(ct.c_float),
                                             ct.POINTER(ct.c_float),
                                             ct.POINTER(ct.c_float),
                                             ct.POINTER(ct.c_int),
                                             ct.POINTER(ct.c_int),
                                             ct.c_float,
                                             ct.POINTER(ct.c_int64)]
        self.lib.n2_merger3d_write_out.restype=None
        self.lib.n2_merger3d_write_out.argtypes=[ct.c_int,
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_int64),
                                                 ct.POINTER(ct.c_float),
                                                 ct.POINTER(ct.c_float),
                                                 ct.c_float,
                                                 ct.POINTER(ct.c_char)]

    def n2merger(self,x_in,y_in,merg_len):
        x = np.copy(x_in)
        y = np.copy(y_in)
        w = np.zeros(x_in.size,dtype=int)
        size = np.zeros((1),dtype=int)
        size[0] = x_in.size
        merg_len = np.float(merg_len)
        #print x
        #print y
        #print w
        #print size
        self.lib.n2_merger_double(x.ctypes.data_as(ct.POINTER(ct.c_double)),
                                 y.ctypes.data_as(ct.POINTER(ct.c_double)),
                                 w.ctypes.data_as(ct.POINTER(ct.c_int)),
                                 size.ctypes.data_as(ct.POINTER(ct.c_int)),
                                 merg_len)
        #print x
        #print y
        #print w
        #print size
        size_final = size[0]
        return x[:size_final],y[:size_final],w[:size_final]
                                               
    def n2merger_writeout(self,step,merg_len,
                          file_loc,
                          x_in,y_in,z_in,
                          vx_in,vy_in,vz_in,
                          core_tags,
                          infall_mass,
                          infall_time):
        x = np.copy(x_in).astype('f4')
        y = np.copy(y_in).astype('f4')
        z = np.copy(z_in).astype('f4')
        vx = np.copy(vx_in).astype('f4')
        vy = np.copy(vy_in).astype('f4')
        vz = np.copy(vz_in).astype('f4')
        tags = np.copy(core_tags).astype('i8')
        mass = np.copy(infall_mass).astype('f4')
        time = np.copy(infall_time).astype('f4')
        size = int(x.size)

        self.lib.n2_merger3d_write_out(size,
                                       x.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       y.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       z.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       vx.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       vy.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       vz.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       tags.ctypes.data_as(ct.POINTER(ct.c_int64)),
                                       time.ctypes.data_as(ct.POINTER(ct.c_float)),                                     
                                       mass.ctypes.data_as(ct.POINTER(ct.c_float)),
                                       merg_len,
                                       "%s/p_core%d.txt"%(file_loc,step))
                                   
    def n2merger3d(self,x_in,y_in,z_in,merg_len,colors_out=False):
        x = np.copy(x_in).astype('f4')
        y = np.copy(y_in).astype('f4')
        z = np.copy(z_in).astype('f4')
        w = np.zeros(x_in.size,dtype='i4')
        size = np.zeros((1),dtype='i4')
        colors = np.zeros(x_in.size,dtype='i8')

        
        size[0] = x_in.size
        merg_len = np.float(merg_len)
        #self.lib.n2_merger_double3d(x.ctypes.data_as(ct.POINTER(ct.c_double)),
        #                            y.ctypes.data_as(ct.POINTER(ct.c_double)),
        #                            z.ctypes.data_as(ct.POINTER(ct.c_double)),
        #                            w.ctypes.data_as(ct.POINTER(ct.c_int32)),
        #                            size.ctypes.data_as(ct.POINTER(ct.c_int32)),
        #                            merg_len)
        self.lib.n2_merger_float3d(x.ctypes.data_as(ct.POINTER(ct.c_float)),
                                   y.ctypes.data_as(ct.POINTER(ct.c_float)),
                                   z.ctypes.data_as(ct.POINTER(ct.c_float)),
                                   w.ctypes.data_as(ct.POINTER(ct.c_int32)),
                                   size.ctypes.data_as(ct.POINTER(ct.c_int32)),
                                   merg_len,
                                   colors.ctypes.data_as(ct.POINTER(ct.c_int64)))

        size_final = size[0]
        x2 = x[:size_final]
        y2 = y[:size_final]
        z2 = z[:size_final]
        w2 = w[:size_final]
        if(colors_out==False):
            return x2,y2,z2,w2
        else:
            return x2,y2,z2,w2,colors

class ZMRIndexValidator:
    def __init__(self,zmr):
        self.z_bins = zmr['z_bins']
        self.m_bins = zmr['mass_bins']
        self.r_bins = zmr['rad_bins']
        self.z_max = max(self.z_bins)
        self.z_min = min(self.z_bins)
        self.m_max = max(self.m_bins)
        self.m_min = min(self.m_bins)
        self.r_max = max(self.r_bins)
        self.r_min = min(self.r_bins)
        self.rad_bin_area = np.pi*(self.r_bins[1:]**2-self.r_bins[:-1]**2)
        self.r_bin_avg = (self.r_bins[:-1]+self.r_bins[1:])/2.0


    def legal_cluster(self,z,mass):
        if( (z_min<z) and (z<z_max) and
            (m_min<mass) and (m<mass_max)):
            return True
        else:
            return False

    def legal_r(self,rad):
        if( (r_min<rad) and (rad<r_max)):
            return True
        else:
            return False
    
    def get_zi_mi(self,z,mass):
        z = np.atleast_1d(np.array(z))
        mass = np.atleast_1d(np.array(mass))
        z_i = np.digitize(z,self.z_bins)-1
        m_i = np.digitize(mass,self.m_bins)-1
        return z_i[0],m_i[0]

    def get_ri(self,r):
        #r = np.array(r).atleast_1d()
        r_i = np.digitize(r,self.r_bins)-1
        return r_i
    def get_r_bin(self,r):
        H,_ = np.histogram(r,bins=self.r_bins)
        return H

    def get_rad_bin_area(self):
        return self.rad_bin_area

    def make_empty_zmr(self):
        new_zmr = {}
        new_zmr['z_bins']= np.copy(self.z_bins)
        new_zmr['rad_bins']= np.copy(self.r_bins)
        new_zmr['mass_bins']= np.copy(self.m_bins)
        new_zmr['zmr_gal_cnt']  = np.zeros((new_zmr['z_bins'].size-1,
                                            new_zmr['mass_bins'].size-1,
                                            new_zmr['rad_bins'].size-1))
        new_zmr['zmr_gal_cnt_err'] = np.zeros_like(new_zmr['zmr_gal_cnt'])
        new_zmr['zmr_gal_density'] = np.zeros_like(new_zmr['zmr_gal_cnt'])
        new_zmr['zmr_gal_density_err'] = np.zeros_like(new_zmr['zmr_gal_cnt'])
        new_zmr['zmr_gal_density_var'] = np.zeros_like(new_zmr['zmr_gal_cnt'])
        new_zmr['zm_ngal'] = np.zeros((new_zmr['z_bins'].size-1,
                                         new_zmr['mass_bins'].size-1))
        new_zmr['zm_ngal_err'] = np.zeros((new_zmr['z_bins'].size-1,
                                         new_zmr['mass_bins'].size-1))
        new_zmr['zm_counts'] = np.zeros((new_zmr['z_bins'].size-1,
                                         new_zmr['mass_bins'].size-1))

        return new_zmr


class Cluster:
    def __init__(self,z,m200,r200,x0,y0,z0,zmr_validator):
        self.z = z
        self.a = 1.0/(z+1.0)
        self.m200 = m200
        self.r200_phys = r200
        self.r200_comv = r200/self.a
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.core_x = None
        self.core_y = None
        self.core_z = None
        self.core_r = None
        self.core_m = None
        self.zmr_val = zmr_validator

    def set_cores(self,x,y,z,rad,mass,min_pos = 0,max_pos=256.0,ignore_periodic=False):
        if(x.size == 0):
            return
        x_diff = x-self.x0
        y_diff = y-self.y0
        z_diff = z-self.z0
        if(not ignore_periodic):
            for i in range(0,x.size):
                if(x_diff[i] > max_pos/2.0):
                    x[i] -= max_pos
                elif(x_diff[i] < -max_pos/2.0):
                    x[i] += max_pos
                if(y_diff[i] > max_pos/2.0):
                    y[i] -= max_pos
                elif(y_diff[i] < -max_pos/2.0):
                    y[i] += max_pos
                if(z_diff[i] > max_pos/2.0):
                    z[i] -= max_pos
                elif(z_diff[i] < -max_pos/2.0):
                    z[i] += max_pos
        self.core_x = np.atleast_1d(x)
        self.core_y = np.atleast_1d(y)
        self.core_z = np.atleast_1d(z)
        self.core_r = np.copy( rad)
        self.core_m = np.copy(mass)
        self.cores_pos = np.zeros((x.size,3))
        self.cores_pos[:,0] = np.array(x)
        self.cores_pos[:,1] = np.array(y)
        self.cores_pos[:,2] = np.array(z)
        self.cores_xy = np.zeros((x.size,2))
        self.cores_xy[:,0] = np.array(x)
        self.cores_xy[:,1] = np.array(y)
        self.kdtree = scipy.spatial.KDTree(self.cores_pos)
        self.kdtree_xy = scipy.spatial.KDTree(self.cores_xy,3)
        if(self.m200 > 1e15):
            indx = np.argmax(self.core_m)
            print "core_m/r", self.core_m[indx], self.core_r[indx]
        return
       
    def compute_dist_mat(self,x,y):
        dist = np.zeros((x.size,x.size),dtype='float')
        for i in range(0,x.size):
            for j in range(0,i):
                dist[i,j] = (x[i]-x[j])**2+(y[i]-y[j])**2
        return dist

    def merge_colors(self,color,c1,c2):
        color[color==c2]=c1
        return
    
    def group_by_id(self,data):
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

    def merge_cores(self,merg_len,core_x,core_y):
        if( core_x==None or core_x.size==0):
            return np.zeros(0),np.zeros(0),np.zeros(0,dtype='i4')
        if(merg_len==0):
            return core_x,core_y,np.ones_like(core_x,dtype='i4')
        color = np.arange(0,core_x.size)
        merg_len = merg_len**2
        dist_xy = self.compute_dist_mat(core_x,core_y)
        for i in range(0,core_x.size):
            for j in range(0,i):
                dist = (core_x[i]-core_x[j])**2+(core_y[i]-core_y[j])**2
                #dist = self.dist_xx[i,j]+self.dist_yy[i,j]
                #dist = self.dist_xy[i,j]
                if(dist <= merg_len and color[i] != color[j]):
                    self.merge_colors(color,color[i],color[j])
        indx = np.argsort(color)
        color = color[indx]
        x = core_x[indx]
        y = core_y[indx]
        [group_start,group_size] = self.group_by_id(color)
        res_x = np.zeros(group_start.size)
        res_y = np.zeros(group_start.size)
        res_weight = np.zeros(group_start.size)
        for i in range(0,group_start.size):
            res_x[i] = np.average(x[group_start[i]:group_start[i]+group_size[i]])
            res_y[i] = np.average(y[group_start[i]:group_start[i]+group_size[i]])
            res_weight[i] = group_size[i]
        return res_x,res_y,res_weight

    def make_rad_prof(self,dis_len,merg_len):
        res_x, res_y, res_weight = process_cores(dis_len,merg_len)
        
    def process_cores(self,dis_len,merg_len):
        slct = self.core_r<dis_len
        res_x,res_y,res_weight = self.merge_cores(merg_len,
                                                  self.core_x[slct],
                                                  self.core_y[slct])
        return res_x,res_y,res_weight

    def set_n2merg(self,lib):
        self.n2merg = lib

    def process_cores_lib(self,dis_len,merg_len,infall_mass,ndim_3d=False,writeout=False):
        if(infall_mass != None):
            #only disrupted cores 
            slct1 = self.core_r>dis_len
            slct2 = self.core_m>10.0**infall_mass
            slct = slct1 & slct2
            if(ndim_3d):
                return self.n2merg.n2merger3d(self.core_x[slct],self.core_y[slct],self.core_z[slct],merg_len)
            else:
                return self.n2merg.n2merger(self.core_x[slct],self.core_y[slct],merg_len)
        else:
            #only disrupted cores
            slct = self.core_r>dis_len
            if(not ndim_3d):
                return self.n2merg.n2merger(self.core_x[slct],self.core_y[slct],merg_len)
            else:
                return self.n2merg.n2merger3d(self.core_x[slct],self.core_y[slct],self.core_z[slct],merg_len)
    
    def writeout_cores_lib(self,dis_len,merg_len,infall_mass):
        slct1 = self.core_r<dis_len
        slct2 = self.core_m>10.0**infall_mass
        slct = slct1 & slct2
        
    def get_rad_bin_cnt(self,dis_len,merg_len,infall_mass=None,return_pos = False):
        #res_x,res_y,res_weight = self.process_cores_lib(dis_len,merg_len,infall_mass)
        #get BCG like galaxies
        bcg_x,bcg_y,bcg_z,bcg_weight = self.process_cores_lib(dis_len,merg_len,infall_mass,ndim_3d=True)
        bcg_slct = bcg_weight > 1
        bcg_rad = np.sqrt((bcg_x[bcg_slct]-self.x0)**2+(bcg_y[bcg_slct]-self.y0)**2)/self.r200_comv
        #compact objects
        if(infall_mass==None):
            gal_slct = (self.core_r<dis_len) 
            gal_slct2 = (self.core_r>dis_len)
        else:
            gal_slct = (self.core_r<dis_len) & (self.core_m>10**infall_mass)
            gal_slct2 = (self.core_r>dis_len) & (self.core_m>10**infall_mass)

        gal_x = self.core_x[gal_slct];
        gal_y = self.core_y[gal_slct];
        gal_z = self.core_x[gal_slct];
        
        gal_rad = np.sqrt((gal_x-self.x0)**2+(gal_y-self.y0)**2)/self.r200_comv
        r_bin_cnt = self.zmr_val.get_r_bin(bcg_rad)+self.zmr_val.get_r_bin(gal_rad)
        if(self.m200 > 5e16):
            print bcg_x[bcg_slct],bcg_y[bcg_slct]
            plt.figure()
            plt.scatter(gal_x,gal_y,'s',ec='blue',c='None')
            plt.scatter(self.core_x[sgal_slct2],self.core_y[gal_slct2],'s',ec='blue',fc='None')
            plt.plot(bcg_x[bcg_slct],bcg_y[bcg_slct],'go')
            plt.show()
        r_bin_err = np.sqrt(r_bin_cnt)
        ngal = np.sum(gal_rad<self.r200_comv)+np.sum(bcg_rad<self.r200_comv)
        #print "cluster ngal:",ngal
        for i in range(0,len(r_bin_err)):
            if(r_bin_err[i]==0):
                r_bin_err[i] = 1.0
        if(return_pos):
            return r_bin_cnt,r_bin_err,ngal,bcg_x,bcg_y,bcg_z,bcg_weight
        else:
            return r_bin_cnt,r_bin_err,ngal

    def get_rad_prof(self,dis_len,merg_len):
        r_bin_cnt,r_bin_err = self.get_rad_bins(dis_len,merg_len)
        r_bin_area = zmr_val.get_rad_bin_area
        return r_bin_cnt/r_bin_area,r_bin_err/r_bin_area
            
    def get_zi_mi(self):
        return self.zmr_val.get_zi_mi(self.z,self.m200)

    def process_cores_kdtree(self,dis_len,merg_len):
        slct = self.core_r<dis_len
        colors = np.arange(0,len(self.core_x))
        for i in range(0,len(self.core_x)):
            if(slct[i]):
                indx = self.kdtree_xy.query_ball_point((self.core_x[i],self.core_y[i]),merg_len)

        #pairs = self.kdtree_xy.query_pairs(merg_len)


def zmr_from_clusters(dis_len,merg_len,clstrs,zmr_validator,infall_mass=None):
    zmr_core = zmr_validator.make_empty_zmr()
    for clstr in clstrs:
        rad_cnt,rad_cnt_err,ngal = clstr.get_rad_bin_cnt(dis_len,merg_len,infall_mass)
        z_i, m_i = clstr.get_zi_mi()
        zmr_core['zmr_gal_cnt'][z_i,m_i]+=rad_cnt
        zmr_core['zm_counts'][z_i,m_i]+=1.0
        zmr_core['zm_ngal'][z_i,m_i]+=ngal
    zmr_count_to_density(zmr_core,0)
    return zmr_core

num_i = 2
num_j = 3
num_k = 0
def write_param_array(f,name,data):
    f.write(name)
    for i in range(0,len(data)):
        f.write(" "+str(data[i]))
    f.write("\n")
    return

def write_param_value(f,name,data):
    f.write(name+" "+str(data)+"\n")
    return 

def write_out_gal_clusters(clstrs,m_infall,r_disrupt,r_merger):
    i=0
    for clstr in clstrs:
        rad_cnt,rad_cnt_err,ngal,gal_x,gal_y,gal_z,gal_w = clstr.get_rad_bin_cnt(r_disrupt,r_merger,m_infall,return_pos=True)
        f=open('tmp/clstr%d_py.param'%i,'w')
        write_param_value(f,"halo_x",clstr.x0)
        write_param_value(f,"halo_y",clstr.y0)
        write_param_value(f,"halo_z",clstr.z0)
        write_param_value(f,"halo_mass",clstr.m200)
        write_param_value(f,"halo_radius",clstr.r200_comv)
        write_param_value(f,"halo_radius_phys",clstr.r200_phys)
        write_param_value(f,"halo_radius_comv",clstr.r200_comv)
        write_param_array(f,"core_x",clstr.core_x)
        write_param_array(f,"core_y",clstr.core_y)
        write_param_array(f,"core_z",clstr.core_z)
        write_param_array(f,"core_r",clstr.core_r)
        write_param_array(f,"core_m",clstr.core_m)
        write_param_array(f,"gal_x",gal_x)
        write_param_array(f,"gal_y",gal_y)
        write_param_array(f,"gal_z",gal_z)
        write_param_array(f,"gal_w",gal_w)
        write_param_array(f,"gal_type",np.zeros_like(gal_w))
        write_param_array(f,"cprtcl_x",[])
        write_param_array(f,"cprtcl_y",[])
        write_param_array(f,"cprtcl_z",[])
        write_param_array(f,"dis_cp_x",[])
        write_param_array(f,"dis_cp_y",[])
        write_param_array(f,"dis_cp_z",[])
        write_param_array(f,"r_bins",clstr.zmr_val.r_bins)
        write_param_array(f,"r_cnt",rad_cnt)
        write_param_value(f,"Ngal",ngal)
        write_param_value(f,"m_infall",10**m_infall)
        write_param_value(f,"r_disrupt",r_disrupt)
        write_param_value(f,"r_merger",r_merger)
        i+=1

def calc_gal_density_cost(core_zmr,sdss_zmr):
    result =0.0
    z_size = sdss_zmr['z_binns'].size -1
    m_size = sdss_zmr['mass_bins'].size -1
    r_size = sdss_zmr['rad_bins'].size -1
    for i in range(0,z_size):
        for j in range(0,m_size):
            for k in range(0,r_size):
                diff = np.abs(core_zmr['zmr_gal_density'][i,j,k] - sdss_zmr['rad_prof'][i,j,k])
                err  = core_zmr['zmr_gal_density_err'][i,j,k]**2 + sdss_zmr['rad_prof_err'][i,j,k]**2
                cost = diff**2/err

                if(core_zmr['zm_counts'][i,j] == 0 or sdss_zmr['zmr_cnt'][i,j,k] == 0):
                    cost = 0.0


                result += cost
    return result

def calc_gal_density_cost2(core_zmr,sdss_zmr):
    diff = np.square(core_zmr['zmr_gal_density']-sdss_zmr['rad_prof'])
    err = np.square(core_zmr['zmr_gal_density_err'])+np.square(sdss_zmr['rad_prof_err'])
    cost = diff/err
    
    zeros = np.zeros_like(cost)
    cost = np.where(sdss_zmr['zmr_cnt']>=0,cost,zeros)
    z_size = sdss_zmr['z_bins'].size -1
    m_size = sdss_zmr['mass_bins'].size -1
    r_size = sdss_zmr['rad_bins'].size -1
    for i in range(0,z_size):
        for j in range(0,m_size):
            if(core_zmr['zm_counts'][i,j] == 0 or np.sum(cost[i,j,:])>1e20):
                for k in range(0,r_size):
                    cost[i,j,k]=0

    if(np.sum(cost>1e20)):
        print "void this!!"
        for i in range(0,z_size):
            for j in range(0,m_size):
                if(np.sum(cost[i,j,:]>1e20)):
                   print i,j
                   print core_zmr['zmr_gal_density'][i,j,:]
                   print sdss_zmr['rad_prof'][i,j,:]
                   print core_zmr['zmr_gal_density'][i,j,:]
                   print sdss_zmr['rad_prof_err'][i,j,:]
                   print diff[i,j,:]
                   print err[i,j,:]
                   print cost[i,j,:]

        exit()
    return np.sum(cost)

def zmr_diff(param):
    disruption_len = param[0]
    merger_len  = param[1]
    ##print "\n\ntrying dis: ",disruption_len," merger_len: ",merger_len
    core_zmr = make_core_zmr(disruption_len,merger_len)
    cost_mat_gal_den = calc_gal_density_cost(core_zmr,zmr_sdss)
    ##print "cost: ", cost_mat_gal_den
    return cost_mat_gal_den

def zmr_count_to_density(zmr,minimum=0.1):
    rad_bins = zmr['rad_bins']
    rad_bins_avg = (rad_bins[:-1]+rad_bins[1:])/2.0
    rad_bins_area = np.pi*(rad_bins[1:]**2-rad_bins[:-1]**2)
    #making ngal & ngal_err
    #print zmr['zm_ngal']
    #print zmr['zm_counts']
    zmr['zm_ngal_err']=np.sqrt(zmr['zm_ngal'])/zmr['zm_counts']
    zmr['zm_ngal']/=zmr['zm_counts']
    #print zmr['zm_ngal']

    for i in range(0,zmr['z_bins'].size-1):
        for j in range(0,zmr['mass_bins'].size-1):
            #making zmr gal density
            zmr['zmr_gal_cnt_err'][i,j] = np.sqrt(zmr['zmr_gal_cnt'][i,j])
            for k in range(0,zmr['rad_bins'].size-1):
                if(zmr['zmr_gal_cnt_err'][i,j,k] == 0):
                    zmr['zmr_gal_cnt_err'][i,j,k]=1.0
            if(zmr['zm_counts'][i,j]>0):
                zmr['zmr_gal_density'][i,j]=zmr['zmr_gal_cnt'][i,j]/zmr['zm_counts'][i,j]/rad_bins_area
                zmr['zmr_gal_density_err'][i,j]=zmr['zmr_gal_cnt_err'][i,j]/zmr['zm_counts'][i,j]/rad_bins_area

def npzfile_to_dic(npzfile):
    result = {}
    for key in npzfile.keys():
        result[key] = npzfile[key]
    return result

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

def save_processed_core_cat(loc,core_cat,intact_slct,colors):
    dtk.ensure_dir(loc)
    hfile = h5py.File(loc,mode='a')
    steps = core_cat.get_steps()
    print "Making hdf5 file"
    i = 0
    for step in steps:
        print "\tworking on step",step
        step_group = hfile.require_group('%d'%step)
        unique_colors = np.unique(colors)
        core_num = unique_colors.size
        gal_x = np.zeros(core_num,dtype='f4')
        gal_y = np.zeros(core_num,dtype='f4')
        gal_z = np.zeros(core_num,dtype='f4')
        gal_vx = np.zeros(core_num,dtype='f4')
        gal_vy = np.zeros(core_num,dtype='f4')
        gal_vz = np.zeros(core_num,dtype='f4')
        gal_vx2 = np.zeros(core_num,dtype='f4')
        gal_vy2 = np.zeros(core_num,dtype='f4')
        gal_vz2 = np.zeros(core_num,dtype='f4')
        gal_tag = np.zeros(core_num,dtype='i8')
        gal_merg_num = np.zeros(core_num,dtype='i4')
        gal_infall_mass_sum = np.zeros(core_num,dtype='f4')
        gal_infall_mass_min = np.zeros(core_num,dtype='f4')
        gal_infall_mass_max = np.zeros(core_num,dtype='f4')
        gal_infall_time_avg = np.zeros(core_num,dtype='f4')
        gal_infall_time_min = np.zeros(core_num,dtype='f4')
        gal_infall_time_max = np.zeros(core_num,dtype='f4')
        gal_host_id = np.zeros(core_num,dtype='i8')
        x = core_cat[step]['x'][intact_slct]
        y = core_cat[step]['y'][intact_slct]
        z = core_cat[step]['z'][intact_slct]
        vx = core_cat[step]['vx'][intact_slct]
        vy = core_cat[step]['vy'][intact_slct]
        vz = core_cat[step]['vz'][intact_slct]
        core_tag = core_cat[step]['core_tag'][intact_slct]
        infall_mass = core_cat[step]['infall_mass'][intact_slct]
        infall_step = core_cat[step]['infall_step'][intact_slct]
        for i in range(0,unique_colors.size):
            #if(i%100==0):
            print i,"/",unique_colors.size
            c = unique_colors[i]
            slct= colors==c
            #gal_xyz
            gal_x[i] = np.average(x[slct])
            gal_y[i] = np.average(y[slct])
            gal_z[i] = np.average(z[slct])
            #gal_vel_xyz
            gal_vx[i] = np.average(vx[slct])
            gal_vy[i] = np.average(vy[slct])
            gal_vz[i] = np.average(vz[slct])
            #gal_vel_xyz median
            gal_vx2[i] = np.median(vx[slct])
            gal_vy2[i] = np.median(vy[slct])
            gal_vz2[i] = np.median(vz[slct])

            #gal_tag
            gal_tag[i] = np.min(core_tag[slct])
            #gal_merg_num
            gal_merg_num[i] = np.sum(slct)
            #infall mass
            infall_masses = infall_mass[slct]
            gal_infall_mass_sum[i] = np.sum(infall_masses)
            gal_infall_mass_min[i] = np.min(infall_masses)
            gal_infall_mass_max[i] = np.max(infall_masses)
            #infall times
            infall_times = infall_step[slct]
            gal_infall_time_avg[i] = np.average(infall_times)
            gal_infall_time_min[i] = np.min(infall_times)
            gal_infall_time_max[i] = np.max(infall_times)
            #gal_host_id
            gal_host_id[i] = stats.mode(core_cat[step]['fof_halo_tag'])[0]
            
            #putting the above data into the hdf5 file
        step_group.create_dataset('x',data=gal_x)
        step_group.create_dataset('y',data=gal_y)
        step_group.create_dataset('z',data=gal_z)
        step_group.create_dataset('vx',data=gal_vx)
        step_group.create_dataset('vy',data=gal_vy)
        step_group.create_dataset('vz',data=gal_vz)
        step_group.create_dataset('vx_median',data=gal_vx2)
        step_group.create_dataset('vy_median',data=gal_vy2)
        step_group.create_dataset('vz_median',data=gal_vz2)
        step_group.create_dataset('core_tag',data=gal_tag)
        step_group.create_dataset('infall_mass',data=gal_infall_mass_sum)
        step_group.create_dataset('infall_mass_min',data=gal_infall_mass_min)
        step_group.create_dataset('infall_mass_max',data=gal_infall_mass_max)
        step_group.create_dataset('infall_time',data=gal_infall_time_avg)
        step_group.create_dataset('infall_time_min',data=gal_infall_time_min)
        step_group.create_dataset('infall_time_max',data=gal_infall_time_max)
    hfile.close()
    return







