#!/usr/bin/env python

#Author: Duncan Campbell
#March, 2015
#Yale University
#examine the quenched fraction of neighbors as a function of seperation

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from halotools.mock_observables import clustering

def main():

    sm_min=9.5
    sm_max=10.0

    #open mock
    #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
    catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle'
    #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_sat_shuffle'
    #catalogue = 'sm_8.5_s0.2_sfr_c-1.0_125'
    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    Lbox = np.array([250.0, 250.0, 250.0])
    
    #define centrals and satellites
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    
    #define quenched and star-forming
    LHS = -11.0
    blue = (mock['SSFR']>LHS) #indices of blue galaxies
    red = (mock['SSFR']<LHS) #indicies of red galaxies
    
    #define samples
    red_host = host & red
    blue_host = host & blue
    
    #define sm threshold for central galaxies
    sm_bin = (mock['Mstar']>sm_min) & (mock['Mstar']<sm_max)
    red_host = red_host & sm_bin
    blue_host = blue_host & sm_bin
    
    #count pairs
    rbins = np.linspace(-1,1.25,20)
    rbins = 10**rbins
    rbin_centers = (rbins[1:]+rbins[:-1])/2.0
    
    data0 = np.vstack((mock['x'], mock['y'], mock['z'])).T
    data1 = np.vstack((mock['x'], mock['y'], mock['z'])).T[red_host]
    data2 = np.vstack((mock['x'], mock['y'], mock['z'])).T[red]
    data3 = np.vstack((mock['x'], mock['y'], mock['z'])).T[blue]
    
    weights = mock['SSFR']
    
    result = clustering.marked_tpcf(data0, rbins, marks1=weights, wfunc=10, period=Lbox)
    
    print(result)
    
    plt.figure()
    plt.plot(rbin_centers, result)
    plt.xscale('log')
    plt.show()
    
    
if __name__ == '__main__':
    main()