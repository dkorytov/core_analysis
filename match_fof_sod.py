#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk
from catalog_reader import Catalog,frag_to_real


print "test"


fof_cat= Catalog("/home/dkorytov/tmp/mira_uni/m011-499.fofproperties")
sod_cat= Catalog("/home/dkorytov/tmp/mira_uni/m011-499.sodproperties")
#fof_cat= Catalog("/home/dkorytov/tmp/mira_uni2/m011-300.fofproperties")
#sod_cat= Catalog("/home/dkorytov/tmp/mira_uni2/m011-300.sodproperties")
halo_cat = Catalog()
fof_cat.add_step(499)
sod_cat.add_step(499)

fof_cat.add_data_name("fof_halo_tag")
fof_cat.add_data_name("fof_halo_count")
fof_cat.add_data_name("fof_halo_mass")


sod_cat.add_data_name("fof_halo_tag")
sod_cat.add_data_name("fof_halo_count")
sod_cat.add_data_name("sod_halo_mass")
sod_cat.add_data_name("fof_halo_center_x")
sod_cat.add_data_name("fof_halo_center_y")
sod_cat.add_data_name("fof_halo_center_z")

print "reading in the data"
fof_cat.read_gio()
sod_cat.read_gio()


print "merging"
halo_cat.join(fof_cat,sod_cat,join_on='fof_halo_tag',verbose=True)


mass_bins = np.logspace(2,6,20)
mass_bins_avg = (mass_bins[:-1]+mass_bins[1:])/2.0

left_fof,_   = np.histogram(fof_cat[499]['fof_halo_count'],bins=mass_bins)
left_sod,_   = np.histogram(sod_cat[499]['fof_halo_count'],bins=mass_bins)
match_halo,_ = np.histogram(halo_cat[499]['fof_halo_count'],bins=mass_bins)


plt.figure()
plt.title('fof left over')
plt.plot(mass_bins_avg,left_fof,label='fof')
plt.xscale('log')
plt.yscale('log')
plt.plot(mass_bins_avg,left_sod,label='sod')
plt.plot(mass_bins_avg,match_halo,label='both')
plt.legend()

plt.figure()
plt.plot(sod_cat[499]['fof_halo_center_x'],sod_cat[499]['fof_halo_center_y'],'x')

plt.show()
