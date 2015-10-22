#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#make plot of stellar mass vs halo for star-forming and quenched galaxies

#load packages
from __future__ import print_function, division
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import sys
import h5py


def main():

    #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle'
    #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_sat_shuffle'
    #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
    catalogue = 'sm_9.5_s0.2_sfr_c1.0_Halfmass_Scale_250'
        
    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(mock.dtype.names)
    
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    red = (mock['red']==1)
    blue = (mock['red']==0)
    
    halo_prop = 'mvir'
    
    bins = np.arange(9.5,12.5,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    mean_red_halo_mass = cu.statistics.binned_mean(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    mean_blue_halo_mass = cu.statistics.binned_mean(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    std_red_halo_mass = cu.statistics.binned_std(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    std_blue_halo_mass = cu.statistics.binned_std(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    print(std_red_halo_mass*10**mean_red_halo_mass)
    
    #plot stellar mass halo mass relation for the mock
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    p = plt.scatter(10.0**mock['Mstar'][host],mock['mvir'][host],c=mock['SSFR'][host],\
                    cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',
                    lw=0, s=2, rasterized=True)
    #p = plt.plot(10.0**mock['Mstar'][host & red],mock['mvir'][host & red], '.', color='red', ms=2)
    #p = plt.plot(10.0**mock['Mstar'][host & blue],mock['mvir'][host & blue], '.', color='blue', ms=2)
    p1a = plt.errorbar(10**bin_centers, 10**mean_red_halo_mass, yerr=std_red_halo_mass*10**mean_red_halo_mass, color='red')
    p1b = plt.errorbar(10**bin_centers, 10**mean_blue_halo_mass, yerr=std_blue_halo_mass*10**mean_blue_halo_mass, color='blue')
    plt.title('age-matching')
    plt.ylabel(r'$M_{\rm halo}~[h^{-1}M_{\odot}]$')
    plt.xlabel(r'$M_{*} ~[h^{-2}M_{\odot}]$')
    plt.xlim([10.0**9.5,10.0**12.5])
    plt.ylim([10.0**11.0,10.0**15.0])
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
    fig.savefig('/Users/duncan/Desktop/SMHM_age_matching.pdf')

if __name__ == '__main__':
    main()