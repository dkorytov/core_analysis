
from __future__ import print_function


import numpy as np
import dtk
import h5py

class ZMR:
    def load_param(self, file_loc):
        print(file_loc)
        if ".param" in file_loc:
            param = dtk.Param(file_loc)
            self.z_bins = np.array(param.get_float_list("z_bins"))
            self.m_bins = np.array(param.get_float_list("m_bins"))
            self.r_bins = np.array(param.get_float_list("r_bins"))
            self.z_size = self.z_bins.size -1
            self.m_size = self.m_bins.size -1
            self.r_size = self.r_bins.size -1
            zm_shape = (self.z_size,self.m_size)
            zmr_shape = (self.z_size,self.m_size,self.r_size)
            self.zm_Ngal = np.array(param.get_float_list("zm_Ngal")).reshape(zm_shape)
            self.zm_Ngal_err = np.array(param.get_float_list("zm_Ngal_err")).reshape(zm_shape)
            self.zm_Ngal_var = np.array(param.get_float_list("zm_Ngal_var")).reshape(zm_shape)
            self.zmr_gal_counts = np.array(param.get_float_list("zmr_gal_counts")).reshape(zmr_shape)
            self.zmr_gal_density = np.array(param.get_float_list("zmr_gal_density")).reshape(zmr_shape)
            self.zmr_gal_density_err = np.array(param.get_float_list("zmr_gal_density_err")).reshape(zmr_shape)
            self.zmr_gal_density_var = np.array(param.get_float_list("zmr_gal_density_var")).reshape(zmr_shape)
            self.zmr_dgal_dr = np.array(param.get_float_list("zmr_dgal_dr")).reshape(zmr_shape)
            self.zmr_dgal_dr_err = np.array(param.get_float_list("zmr_dgal_dr_err")).reshape(zmr_shape)
            self.zmr_dgal_dr_var = np.array(param.get_float_list("zmr_dgal_dr_var")).reshape(zmr_shape)
            self.zmr_gal_accum = np.array(param.get_float_list("zmr_gal_accum")).reshape(zmr_shape)
            self.zmr_gal_accum_err = np.array(param.get_float_list("zmr_gal_accum_err")).reshape(zmr_shape)
            self.zmr_gal_accum_var = np.array(param.get_float_list("zmr_gal_accum_var")).reshape(zmr_shape)

            self.zm_counts = np.array(param.get_float_list("zm_counts")).reshape(zm_shape)
            self.zmr_counts = np.array(param.get_float_list("zmr_counts")).reshape(zmr_shape) 
        elif ".hdf5" in file_loc:
            hfile = h5py.File(file_loc, 'r')
            print(hfile.keys())
            self.z_bins = hfile['z_bins'].value
            self.m_bins = hfile['m_bins'].value
            self.r_bins = hfile['r_bins'].value
            self.z_size = self.z_bins.size -1
            self.m_size = self.m_bins.size -1
            self.r_size = self.r_bins.size -1
            self.zm_Ngal = hfile['zm_Ngal']
            self.zm_Ngal_err = hfile['zm_Ngal_err']
            self.zm_Ngal_var = hfile['zm_Ngal_var']
            
            self.zmr_gal_density = hfile['zmr_gal_density']
            self.zmr_gal_density_var = hfile['zmr_gal_density_var']
            self.zmr_gal_density_err = hfile['zmr_gal_density_err']
            self.zm_counts = hfile['zm_counts']
            self.zmr_counts = hfile['zmr_counts']
        else:
            raise ValueError("{} doesn't have .param or .hdf5".format(file_loc))


    def __init__(self,file_loc=None):
        if(file_loc != None):
            self.load_param(file_loc)
    
    def cal_cost(self,zmr2):
        cost = 0.0
        for zi in range(0,z_bins.size-1):
            for mi in range(0,m_bins.size-1):
                if(zm_counts[zi,mi]>0):
                    err = np.square(zmr_gal_density_err[zi,mi,:]) +np.square( zmr2.zmr_gal_density_err[zi,mi,:])
                    diff = np.square(zmr_gal_density[zi,mi,:]- zmr2.zmr_gal_density[zi,mi,:])
                    cost += diff/err
        return cost

    def save_as_hdf5(self, file_loc):
        hfile = h5py.File(file_loc, 'w')
        hfile['z_bins'] = self.z_bins
        hfile['m_bins'] = self.m_bins
        hfile['r_bins'] = self.r_bins
        hfile['z_size'] = self.z_size
        hfile['m_size'] = self.m_size
        hfile['r_size'] = self.r_size
        hfile['zm_Ngal'] = self.zm_Ngal
        hfile['zm_Ngal_err'] = self.zm_Ngal_err
        hfile['zm_Ngal_var'] = self.zm_Ngal_var
        hfile['zmr_gal_counts'] = self.zmr_gal_counts
        hfile['zmr_gal_density'] = self.zmr_gal_density
        hfile['zmr_gal_density_err'] = self.zmr_gal_density_err
        hfile['zmr_gal_density_var'] = self.zmr_gal_density_var
        hfile['zmr_dgal_dr'] = self.zmr_dgal_dr
        hfile['zmr_dgal_dr_err'] = self.zmr_dgal_dr_err
        hfile['zmr_dgal_dr_var'] = self.zmr_dgal_dr_var
        hfile['zmr_gal_accum'] = self.zmr_gal_accum
        hfile['zmr_gal_accum_err'] = self.zmr_gal_accum_err
        hfile['zmr_gal_accum_var'] = self.zmr_gal_accum_var
        hfile['zm_counts'] = self.zm_counts
        hfile['zmr_counts'] = self.zmr_counts
        hfile.close()

    def get_ngal(self):
        zmr_ngal = np.zeros((self.z_size, self.m_size))
        r_bin_area = np.reshape(np.diff(self.r_bins**2*np.pi), (1, 1, len(self.r_bins)-1))
        zm_ngal  = np.sum(self.zmr_gal_density*r_bin_area, axis=2)
        zm_ngal_err  = np.sqrt(np.sum((self.zmr_gal_density_err*r_bin_area)**2, axis=2))

        print(zm_ngal_err.shape)
        print(zm_ngal.shape)
        return zm_ngal, zm_ngal_err
