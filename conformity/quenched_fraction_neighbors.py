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
from halotools.mock_observables.pair_counters import rect_cuboid_pairs

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
    rbins = np.linspace(-1,1,20)
    rbin_centers = (rbins[1:]+rbins[:-1])/2.0
    rbins = 10**rbins
    data1 = np.vstack((mock['x'], mock['y'], mock['z'])).T[red_host]
    data2 = np.vstack((mock['x'], mock['y'], mock['z'])).T[red]
    data3 = np.vstack((mock['x'], mock['y'], mock['z'])).T[blue]
    
    print(len(data1), len(data2), len(data3))
    red_pairs = rect_cuboid_pairs.npairs(data1, data2, rbins, period=Lbox)
    red_pairs = np.diff(red_pairs)
    blue_pairs = rect_cuboid_pairs.npairs(data1, data3, rbins, period=Lbox)
    blue_pairs = np.diff(blue_pairs)
    all_pairs = red_pairs+blue_pairs
    
    red_quenched_frac = red_pairs/all_pairs
    
    data1 = np.vstack((mock['x'], mock['y'], mock['z'])).T[blue_host]
    data2 = np.vstack((mock['x'], mock['y'], mock['z'])).T[red]
    data3 = np.vstack((mock['x'], mock['y'], mock['z'])).T[blue]
    
    print(len(data1), len(data2), len(data3))
    red_pairs = rect_cuboid_pairs.npairs(data1, data2, rbins, period=Lbox)
    red_pairs = np.diff(red_pairs)
    blue_pairs = rect_cuboid_pairs.npairs(data1, data3, rbins, period=Lbox)
    blue_pairs = np.diff(blue_pairs)
    all_pairs = red_pairs+blue_pairs
    
    blue_quenched_frac = red_pairs/all_pairs
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(10**rbin_centers,red_quenched_frac, color='red')
    plt.plot(10**rbin_centers,blue_quenched_frac, color='blue')
    plt.xlabel(r'$r~(h^{-1}\rm Mpc)$')
    plt.ylabel(r'$f_{\rm quenched}$')
    #plt.ylim([0,1.0])
    plt.ylim([0.4,0.6])
    #plt.xscale('log')
    plt.show()
    
    
    
if __name__ == '__main__':
    main()