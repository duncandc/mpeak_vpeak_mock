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

    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
    catalogues = ['sm_9.5_s0.0_sfr_c-1.0_250','sm_9.5_s0.1_sfr_c-1.0_250',\
                  'sm_9.5_s0.2_sfr_c-1.0_250','sm_9.5_s0.3_sfr_c-1.0_250']
    
    sigmas = [0.0,0.1,0.2,0.3]

    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.xlabel(r'$r~(h^{-1}\rm Mpc)$')
    plt.ylabel(r'$f_{\rm quenched}$')
    plt.title(r'$10.0<\log(M_{*, \rm host}/M_{\odot}h^{-2})<10.5$')
    plt.ylim([0,1])
    plt.xlim([0.1,10])
    plt.xscale('log')
    a = [1.0,0.8,0.6,0.4,0.2]
    #red_colors=['#FF0000','#FF2900','#FF5C00','#FF8F00','#FFCC00']
    #blue_colors=['#0000FF','#0040FF','#0080FF','#00BFFF','#00FFFF']
    p1 = [None]*5
    p2 = [None]*5
    for i, catalogue in enumerate(catalogues):
        #open catalogue
        print(i,catalogue)
        f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
        mock = f.get(catalogue)
        Lbox = np.array([250.0,250.0,250.0])
    
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
        
        #host mass selection
        bin = (mock['Mstar']>10.0) & (mock['Mstar']<10.5)
        red_host = red_host & bin
        blue_host = blue_host & bin
    
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
    
        p1[i], = plt.plot(10**rbin_centers,red_quenched_frac, color='red', alpha=a[i], lw=1)
        p2[i], =plt.plot(10**rbin_centers,blue_quenched_frac, color='blue', alpha=a[i], lw=1)
    p1 = tuple(p1)
    p2 = tuple(p2)
    plt.legend(p1, (r"$\sigma=0.0$",r"$\sigma=0.1$",r"$\sigma=0.2$",r"$\sigma=0.3$"),\
               loc=4, fontsize=10, frameon=False, labelspacing=0.01, title='quenched centrals')
    fig.savefig('/Users/duncan/Desktop/conformity.pdf')
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'conformity_rho_effect'
    fig.savefig(savepath+filename+'.pdf')
    
    
if __name__ == '__main__':
    main()