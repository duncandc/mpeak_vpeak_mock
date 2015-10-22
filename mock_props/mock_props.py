#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#print out some basic mock properties

#load packages
from __future__ import print_function, division
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import sys
import h5py


def main():

    catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
    #catalogue = 'sm_9.5_s0.0_sfr_c-1.0_250_cen_shuffle'

    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(mock.dtype.names)

    hosts = (mock['upid']==-1)
    subs = (mock['upid']!=-1)

    N = len(mock)
    Ncen = np.sum(hosts)
    Nsat = np.sum(subs)

    print('N galaxies: ', N)
    print('N centrals: ', Ncen)
    print('N satellites: ', Nsat)
    print('N ejected centrals: ', np.sum(mock['ejected']))
    print('N abandoned satellites: ',np.sum(mock['abandoned']))
    print('N ghosted satellites: ',np.sum(mock['ghost']))
    
    print('satellite fraction: ', Nsat/(1.0*(Nsat+Ncen)))


if __name__ == '__main__':
    main()