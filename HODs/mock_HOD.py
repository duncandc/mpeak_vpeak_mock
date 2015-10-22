#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#make plot of mock HOD

#load packages
from __future__ import print_function, division
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import sys
import h5py


def main():

    catalogue = 'sm_9.49_s0.2_sfr_c0.0_250'
    #catalogue = 'sm_9.49_s0.0_sfr_c-1.0_250_cen_shuffle'

    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(mock.dtype.names)
    
    #open halo catalogue (used to account for empty haloes)
    sys.path.insert(0, '../mocks/')
    from make_mock import get_mock
    HC, Lbox = get_mock('Bolshoi')
    
    #define host and sub haloes
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    HC_host = (HC['upid']==-1)
    HC_sub = (HC['upid']!=-1)
    
    #star forming and quenched
    LHS = -11.0
    blue = (mock['SSFR']>LHS) #indices of blue galaxies
    red = (mock['SSFR']<LHS) #indicies of red galaxies

    #calculate the number of members in a halo
    from halotools.utils import aggregation
    def function(x):
        return len(x)
    mock = aggregation.add_members_property(mock, 'group_ID', 'N', function)
    
    def function(x):
        inds = (x['SSFR']<=-11.0)
        red = np.zeros(len(x))
        red[inds]=1.0
        return red
    mock = aggregation.add_members_property(mock, 'group_ID', 'red', function)
    
    def function(x):
        count = (x['Mstar']>=9.5)
        return np.sum(x['red'][count])
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_red', function)
    
    def function(x):
        count = (x['Mstar']>=9.5)
        return len(x[count])-np.sum(x['red'][count])
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_blue', function)
    
    def function(x):
        count = (x['Mstar']>=9.5)
        sat = (x['upid']!=-1)
        return np.sum(x['red'][count&sat])
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_red_sat', function)
    
    def function(x):
        count = (x['Mstar']>=9.5)
        sat = (x['upid']!=-1)
        return len(x[count&sat])-np.sum(x['red'][count&sat])
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_blue_sat', function)
    
    #bin galaxies by halo mass
    bins = np.arange(10,15,0.1)
    bins = 10.0**bins
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    result = np.digitize(mock['Mvir_host'],bins)
    N_cen = np.histogram(mock['mvir'][host],bins)[0]
    N_haloes = np.histogram(HC['mvir'][HC_host],bins)[0]
    
    N = np.zeros(len(bins)-1)
    N_red = np.zeros(len(bins)-1)
    N_blue= np.zeros(len(bins)-1)
    Ncen = np.zeros(len(bins)-1)
    Nsat = np.zeros(len(bins)-1)
    Ncen_red = np.zeros(len(bins)-1)
    Nsat_red = np.zeros(len(bins)-1)
    Ncen_blue = np.zeros(len(bins)-1)
    Nsat_blue = np.zeros(len(bins)-1)
    for i in range(0,len(bins)-1):
        inds = (result==i+1)
        mock_c = mock[inds]
        
        #empty haloes
        N_empty = np.float(N_haloes[i]-N_cen[i])
        zeros = np.zeros(N_empty)
        ones = np.ones(N_cen[i])
       
        #red galaxies
        N_red[i] = np.mean(np.hstack((mock_c['N_red'],zeros)))
        Nsat_red[i] = np.mean(np.hstack((mock_c['N_red_sat'],zeros)))
        Ncen_red[i] = np.mean(np.hstack((mock_c['N_red']-mock_c['N_red_sat'],zeros)))
        
        #blue galaxies
        N_blue[i] = np.mean(np.hstack((mock_c['N_blue'],zeros)))
        Nsat_blue[i] = np.mean(np.hstack((mock_c['N_blue_sat'],zeros)))
        Ncen_blue[i] = np.mean(np.hstack((mock_c['N_blue']-mock_c['N_blue_sat'],zeros)))
        
        #all galaxies
        N[i] = np.mean(np.hstack((mock_c['N'],zeros)))
        Ncen[i] = np.mean(np.hstack((ones,zeros)))
        Nsat[i] = np.mean(np.hstack((N[i]-ones,zeros)))
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.plot(bin_centers, N, color='black')
    plt.plot(bin_centers, Ncen, '--',color='black')
    plt.plot(bin_centers, Nsat, ':',color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0.1,100])
    plt.xlim([10**10.5,10**15])
    plt.xlabel(r'$M_{\rm vir} ~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$\langle N \rangle$')
    plt.show()
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.plot(bin_centers, N, color='black')
    plt.plot(bin_centers, N_red, color='red')
    plt.plot(bin_centers, N_blue, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0.1,100])
    plt.xlim([10**10.5,10**15])
    plt.xlabel(r'$M_{\rm vir} ~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$\langle N \rangle$')
    plt.show()
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.plot(bin_centers, N_red, color='red')
    plt.plot(bin_centers, Ncen_red, '--', color='red')
    plt.plot(bin_centers, Nsat_red, ':', color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0.1,100])
    plt.xlim([10**10.5,10**15])
    plt.xlabel(r'$M_{\rm vir} ~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$\langle N \rangle$')
    plt.show()
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.plot(bin_centers, N_blue, color='blue')
    plt.plot(bin_centers, Ncen_blue, '--', color='blue')
    plt.plot(bin_centers, Nsat_blue, ':', color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([0.1,100])
    plt.xlim([10**10.5,10**15])
    plt.xlabel(r'$M_{\rm vir} ~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$\langle N \rangle$')
    plt.show()

if __name__ == '__main__':
    main()