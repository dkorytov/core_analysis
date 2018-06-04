#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import dtk


print 'reading'
step = 487
#tags = dtk.gio_read('/media01_12_17.AlphaQ.'+str(step)+'.treenodes','fof_halo_tag')
#mass = dtk.gio_read('/home/dkorytov/tmp/01_12_17.AlphaQ.'+str(step)+'.treenodes','fof_halo_mass')
tags = dtk.gio_read('/media/luna1/rangel/AlphaQ/temp/01_12_17.AlphaQ.'+str(step)+'.treenodes','fof_halo_tag')
mass = dtk.gio_read('/media/luna1/rangel/AlphaQ/temp/01_12_17.AlphaQ.'+str(step)+'.treenodes','fof_halo_mass')

print tags
print np.min(tags),np.max(tags)
print 'done reading'
mass_bins = np.logspace(10,16,20)
mass_bins_avg =(mass_bins[:-1]+mass_bins[1:])/2.0
lim = 1024*1024*1024
fof_slct = tags<lim
frag_slct = tags>lim
print lim
print np.sum(fof_slct )
print np.sum(frag_slct)
fof,_ = np.histogram(mass[fof_slct],bins=mass_bins)
frag,_ = np.histogram(mass[frag_slct],bins=mass_bins)
tot = frag+fof
print 'done histing'

plt.figure()
plt.plot(mass_bins_avg,fof.astype('f4')/tot.astype('f4'))
plt.show()
