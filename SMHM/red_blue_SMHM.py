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

    if len(sys.argv)==2:
        rho = 'NA'
        sigma = 'NA'
        catalogue = sys.argv[1]
    elif len(sys.argv)>2:
        rho = sys.argv[1]
        sigma = sys.argv[2]
        catalogue = 'sm_9.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_250'
    else:
        #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle'
        #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_sat_shuffle'
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        rho = -1.0
        sigma = 0.2
        
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
    
    bins = np.arange(9.5,12.5,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    mean_red_halo_mass = cu.statistics.binned_mean(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    mean_blue_halo_mass = cu.statistics.binned_mean(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    std_red_halo_mass = cu.statistics.binned_std(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    std_blue_halo_mass = cu.statistics.binned_std(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    
    #more 2011 results
    red_M0 = 10**10.84
    red_M1 = 10**12.18
    red_gamma_1 = 3.34
    red_gamma_2 = 0.22
    blue_M0 = 10**9.38
    blue_M1 = 10**11.32
    blue_gamma_1 = 2.41
    blue_gamma_2 = 1.12
    
    msample = np.logspace(10,15,100)
    red_M = red_M0 *((msample/red_M1)**red_gamma_1/(1+(msample/red_M1))**(red_gamma_1-red_gamma_2))
    blue_M = blue_M0 *((msample/blue_M1)**blue_gamma_1/(1+(msample/blue_M1))**(blue_gamma_1-blue_gamma_2))
    
    #plot stellar mass halo mass relation for the mock
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.78, bottom=0.2, top=0.9)
    plt.plot(10**bin_centers,10**mean_red_halo_mass,color='red')
    plt.plot(10**bin_centers,10**mean_blue_halo_mass,color='blue')
    plt.fill_between(10**bin_centers, 10**(mean_red_halo_mass+std_red_halo_mass),\
                                      10**(mean_red_halo_mass-std_red_halo_mass),\
                                      color='red', alpha=0.5)
    plt.fill_between(10**bin_centers, 10**(mean_blue_halo_mass+std_blue_halo_mass),\
                                      10**(mean_blue_halo_mass-std_blue_halo_mass),\
                                      color='blue', alpha=0.5)
    plt.plot(red_M, msample, '--', color='red')
    plt.plot(blue_M, msample, '--', color='blue')
    plt.text(10.0**9.6,10**13.8,r'$\rho_{\rm SSFR}=$'+str(rho))
    plt.text(10.0**9.6,10**13.7,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.ylabel(r'$M_{\rm halo}~[h^{-1}M_{\odot}]$')
    plt.xlabel(r'$M_{*} ~[h^{-2}M_{\odot}]$')
    plt.xlim([10.0**9.5,10.0**12.5])
    plt.ylim([10.0**11.5,10.0**14.0])
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

if __name__ == '__main__':
    main()