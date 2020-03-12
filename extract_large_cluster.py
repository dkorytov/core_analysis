#!/usr/bin/env python3



import numpy as np
import matplotlib.pyplot as plt
import h5py


def extract_large_cluster(input_fname='tmp_hdf5/clusters_OR_M200c.hdf5', output_fname='tmp_hdf5/single_cluster.hdf5'):
    if input_fname == output_fname:
        raise
    hfile = h5py.File(input_fname,'r')
    cluster_masses = hfile['cluster/sod_mass'][()]

    print(cluster_masses)
    index = np.argmax(cluster_masses)
    cluster_x = hfile['cluster/x'][index]
    cluster_y = hfile['cluster/y'][index]
    cluster_z = hfile['cluster/z'][index]
    print(np.log10(cluster_masses[index]))
    offset = hfile['cluster/core_offset'][index]
    size = hfile['cluster/core_size'][index]
    core_x = hfile['cores/core_x'][offset:offset+size]
    core_y = hfile['cores/core_y'][offset:offset+size]
    core_z = hfile['cores/core_z'][offset:offset+size]
    core_m = hfile['cores/core_m'][offset:offset+size]
    plt.figure()
    plt.plot(core_x, core_y, '.')
    
    plt.figure()
    plt.plot(core_x, core_z, '.')
    # sort it by infall mass
    srt = np.argsort(-core_m)
    
    hfile_out = h5py.File(output_fname,'w')
    hfile_out['x'] = core_x[srt]
    hfile_out['y'] = core_y[srt]
    hfile_out['z'] = core_z[srt]

    hfile_out.close()

if __name__ == "__main__":
    extract_large_cluster()
    plt.show()
