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
    
    bins = np.arange(9.5,12.5,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    mean_red_halo_mass_age_matching = cu.statistics.binned_mean(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    mean_blue_halo_mass_age_matching = cu.statistics.binned_mean(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    std_red_halo_mass_age_matching = cu.statistics.binned_std(mock[host&red], bins, 'Mstar', halo_prop, use_log=True)[1]
    std_blue_halo_mass_age_matching = cu.statistics.binned_std(mock[host&blue], bins, 'Mstar', halo_prop, use_log=True)[1]
    
    #more 2011 results
    red_M0 = 10**10.84
    red_M1 = 10**12.18
    red_gamma_1 = 3.34
    red_gamma_2 = 0.22
    blue_M0 = 10**9.38
    blue_M1 = 10**11.32
    blue_gamma_1 = 2.41
    blue_gamma_2 = 1.12
    
    msample = np.logspace(10.0,14.0,20)
    red_M = red_M0 *((msample/red_M1)**red_gamma_1/(1+(msample/red_M1))**(red_gamma_1-red_gamma_2))
    blue_M = blue_M0 *((msample/blue_M1)**blue_gamma_1/(1+(msample/blue_M1))**(blue_gamma_1-blue_gamma_2))
    
    #Mandelbaum
    Mstar_blue = np.array([9.93968, 10.283493, 10.589262, 10.835646, 11.004268, 11.191629, 11.338512])
    Mhalo_blue = np.array([11.816899, 11.770057, 12.1783905, 12.594514, 12.703571, 12.821066, 12.814001])
    Mhalo_blue_err = np.array([[11.985507,11.612837],[11.877847,11.57971],[12.242236,12.05176],\
                               [12.68944,12.47412],[12.871635,12.424431],[13.236025,11.803312],[13.385093,11.0]])
    Mstar_red = np.array([10.042392, 10.348698, 10.623265, 10.850929, 11.028584, 11.209333, 11.402584])
    Mhalo_red = np.array([12.1738615, 12.126708, 12.493375, 12.884498, 13.283488, 13.70735, 14.1])
    
    print(Mhalo_blue_err[:,0])
    Mhalo_blue_err = 10.0**Mhalo_blue_err
    Mhalo_blue_err[:,0] = np.fabs(Mhalo_blue_err[:,0] - 10**Mhalo_blue)
    Mhalo_blue_err[:,1] = np.fabs(Mhalo_blue_err[:,1] - 10**Mhalo_blue)
    
    Mhalo_blue_err = np.flipud(Mhalo_blue_err.T)
    print(np.shape(Mhalo_blue_err))
    
    #plot stellar mass halo mass relation for the mock
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    p1a, = plt.plot(10**bin_centers,10**mean_red_halo_mass,color='red')
    p1b, = plt.plot(10**bin_centers,10**mean_blue_halo_mass,color='blue')
    p2a, = plt.plot(10**bin_centers,10**mean_red_halo_mass_age_matching,color='red', linestyle='--', lw=4)
    p2b, = plt.plot(10**bin_centers,10**mean_blue_halo_mass_age_matching,color='blue', linestyle='--', lw=4)
    plt.fill_between(10**bin_centers, 10**(mean_red_halo_mass+std_red_halo_mass),\
                                      10**(mean_red_halo_mass-std_red_halo_mass),\
                                      color='red', alpha=0.5)
    plt.fill_between(10**bin_centers, 10**(mean_blue_halo_mass+std_blue_halo_mass),\
                                      10**(mean_blue_halo_mass-std_blue_halo_mass),\
                                      color='blue', alpha=0.5)
    d1a, = plt.plot(red_M, msample, ':', color='red', lw=4)
    d1b, = plt.plot(blue_M, msample, ':', color='blue', lw=4)
    d2a = plt.errorbar(10**Mstar_blue, 10**Mhalo_blue, yerr=Mhalo_blue_err, fmt='o', color='blue', mec='none')
    d2b, = plt.plot(10**Mstar_red, 10**Mhalo_red, 'o', color='red', mec='none')
    legend2 = plt.legend((d2b,d1a), ('Mandelbaum+ 2015', 'More+ 2011'), loc=2, fontsize=10, frameon=False)
    plt.gca().add_artist(legend2)
    plt.legend((p1a,p2a), ('this work','age-matching'), loc=4, fontsize=10, frameon=False)
    #plt.text(10.0**9.6,10**13.8,r'$\rho_{\rm SSFR}=$'+str(rho))
    #plt.text(10.0**9.6,10**13.7,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.ylabel(r'$M_{\rm halo}~[h^{-1}M_{\odot}]$')
    plt.xlabel(r'$M_{*} ~[h^{-2}M_{\odot}]$')
    plt.xlim([10.0**9.5,10.0**12.5])
    plt.ylim([10.0**11.0,10.0**15.0])
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
    fig.savefig('/Users/duncan/Desktop/SMHM.pdf')

if __name__ == '__main__':
    main()