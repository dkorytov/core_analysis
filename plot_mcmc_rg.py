#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk
import sys





def gelman_ruben_test(data,chain_id):
    ids = np.unique(chain_id)
    n = float(data.size)
    m = float(ids.size)
    ssq = []
    thetab = []
    global_avg = np.mean(data)
    B = 0
    for ch_id in ids:
        slct = chain_id == ch_id
        ssq.append(np.var(data[slct]))
        chain_avg= np.mean(data[slct])
        B += (chain_avg - global_avg)**2
    W = np.mean(ssq)
    B = n/(m-1.0)*B
    var = (n-1.0)/n*W + 1.0/n*B
    R = np.sqrt(var/W)
    return R

if(__name__ == '__main__'):
    print "Gelman Ruben Test"
    param_file_name = sys.argv[1]
    param = dtk.Param(param_file_name)
    file_loc = "output/"+param_file_name+"/mcmc.gio"
    print "loading data"
    mcmc_id = dtk.gio_read(file_loc,"mcmc_walker_id")
    mcmc_step = dtk.gio_read(file_loc,"mcmc_walker_step")
    mcmc_mi = dtk.gio_read(file_loc,"mcmc_mass_infall")
    mcmc_rd = dtk.gio_read(file_loc,"mcmc_r_disrupt")
    
    print "data loaded"
    slct_half = mcmc_step>np.max(mcmc_step)/2.0
    print "M_infall: ", gelman_ruben_test(mcmc_mi[slct_half],mcmc_id[slct_half])
    print "R_dispupt: ", gelman_ruben_test(mcmc_rd[slct_half],mcmc_id[slct_half])
    if(param.get_bool("fit_r_merger")):
        mcmc_rm = dtk.gio_read(file_loc,"mcmc_r_merger")
        print "R_merger: ",gelman_ruben_test(mcmc_mi[slct_half],mcmc_id[slct_half])
    if(param.get_bool("fit_r_fof")):
        mcmc_rm = dtk.gio_read(file_loc,"mcmc_r_fof")
        print "R_fof: ",gelman_ruben_test(mcmc_mi[slct_half],mcmc_id[slct_half])
    
