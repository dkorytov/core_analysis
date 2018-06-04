#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import dtk

params = dtk.Param("mcmc/test.param")


m_i = np.array(params.get_float_list("m_i"))
r_d = np.array(params.get_float_list("r_d"))
r_f = np.array(params.get_float_list("r_f"))
r_m = np.array(params.get_float_list("r_m"))
print min(m_i),max(m_i)
print min(r_d),max(r_d)

plt.figure()
plt.plot(m_i,r_d,'.',alpha=.3)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('m_infall')
plt.ylabel('r_disrupt')
plt.grid();

plt.figure()
H,x_bins,y_bins = np.histogram2d(m_i,r_d,bins=(np.logspace(10,12,100),np.logspace(-2,-1,100)))
print H
plt.pcolor(x_bins,y_bins,H.T,cmap='PuBu',norm = clrs.LogNorm())
plt.yscale('log')
plt.xscale('log')
plt.grid()

plt.figure()
H,x_bins = np.histogram(m_i,bins=(np.logspace(10,12,100)))
xavg = (x_bins[1:]+x_bins[:-1])/2.0
plt.plot(xavg,H,'x-')
plt.xscale('log')
plt.title('m_infall')


plt.figure()
H,x_bins = np.histogram(r_d,bins=(np.logspace(-2,0,100)))
xavg = (x_bins[1:]+x_bins[:-1])/2.0
plt.plot(xavg,H,'x-')
plt.xscale('log')
plt.title("r_dist")

plt.show()

plt.figure()
plt.plot(m_i,r_m,'.',alpha=.3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('m_infall')
plt.ylabel('r_merger')

