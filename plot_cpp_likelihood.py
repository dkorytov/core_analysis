#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk
from matplotlib.colors import LogNorm
import sys

param_file_name = sys.argv[1]
lgrid_param = dtk.Param("output/"+sys.argv[1]+"/lgrid.param")
result = np.array(lgrid_param.get_float_list("result"))
nan_slct = np.isnan(result)
result[nan_slct] = np.ones(np.sum(nan_slct))*1000000
mi_bins = np.array(lgrid_param.get_float_list("mi_bins"))
rd_bins = np.array(lgrid_param.get_float_list("rd_bins"))
rm_bins = np.array(lgrid_param.get_float_list("rm_bins"))

result2 = result.reshape((mi_bins.size,rd_bins.size,rm_bins.size))

lkhd = result2


#mi
mi_lkhd = np.min(lkhd,axis=(1,2))
#rd
rd_lkhd = np.min(lkhd,axis=(0,2))
#rm
rm_lkhd = np.min(lkhd,axis=(0,1))
#mi rd
mi_rd_lkhd = np.min(lkhd,axis=2)
#mi rm
mi_rm_lkhd = np.min(lkhd,axis=1)
#rm rd
rd_rm_lkhd = np.min(lkhd,axis=0)

plt.figure()
print mi_lkhd
plt.plot(mi_bins,mi_lkhd)
plt.xlabel('M_infall')
plt.ylabel('likelihood')
plt.xscale('log')
plt.yscale('log')

plt.figure()
print rd_lkhd

plt.plot(rd_bins,rd_lkhd)
plt.xlabel('R_disrupt')
plt.ylabel('likelihood')
plt.xscale('log')
plt.yscale('log')

plt.figure()
print rm_lkhd
plt.plot(rm_bins,rm_lkhd)
plt.xlabel('R_merger')
plt.ylabel('likelihood')
plt.xscale('log')
plt.yscale('log')

plt.figure()
plt.pcolor(mi_bins,rd_bins,1.0/mi_rd_lkhd.T,norm=LogNorm(),cmap='nipy_spectral')
plt.colorbar()
plt.xscale('log')
plt.xlabel('M_infall')
plt.ylabel('R_disrupt')
plt.xlim((np.min(mi_bins),np.max(mi_bins)))
plt.ylim((np.min(rd_bins),np.max(rd_bins)))

plt.figure()
plt.pcolor(mi_bins,rm_bins,1.0/mi_rm_lkhd.T,norm=LogNorm(),cmap='nipy_spectral')
plt.colorbar()
plt.xscale('log')
plt.xlabel('M_infall')
plt.ylabel('R_merger')

plt.figure()
plt.pcolor(rd_bins,rm_bins,1.0/rd_rm_lkhd.T,norm=LogNorm(),cmap='nipy_spectral')
plt.colorbar()
plt.xlabel('R_disrupt')
plt.ylabel('R_merger')

dtk.save_figs(path='figs/'+param_file_name+"/"+__file__+"/",extension='.png')
plt.show()
