#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk

class ZMR:
    def load_from_param(self,file_loc):
        param = dtk.Param(file_loc)
        self.z_bins = np.array(param.get_float_list("z_bins"))
        self.m_bins = np.array(param.get_float_list("m_bins"))
        self.r_bins = np.array(param.get_float_list("r_bins"))
        self.zmr_gal_density = np.array(param.get_float_list("zmr_gal_density"))
        self.zmr_gal_density_err = np.array(param.get_float_list("zmr_gal_density_err"))
