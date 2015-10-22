#!/usr/bin/env python

#Author: Duncan Campbell
#January 28, 2015
#Yale University
#calculate the satellite quenched fraction in our mocks

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():


    catalogues = ['sm_9.5_s0.0_sfr_c-1.0_250','sm_9.5_s0.1_sfr_c-1.0_250',\
                  'sm_9.5_s0.2_sfr_c-1.0_250','sm_9.5_s0.3_sfr_c-1.0_250']
    sigmas = [0.0,0.1,0.2,0.3]
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    
    alphas = [1.0,0.8,0.6,0.4,0.2]
    p1s =[]
    p2s =[]
    for i, catalogue in enumerate(catalogues): 

        #open mock
        filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
        print('opening mock catalogue:', catalogue+'.hdf5')
        #open catalogue
        f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
        mock = f.get(catalogue)
        mock = np.array(mock)
        print(mock.dtype.names)
    
        #define host haloes and subhaloes
        hosts = (mock['upid']==-1)
        subs = (mock['upid']!=-1)
    
        #define SF and quenched subsamples
        LHS = -11.0
        blue = np.where(mock['SSFR']>LHS)[0]
        red = np.where(mock['SSFR']<LHS)[0]
    
        #calculate satellite fraction
        bins = np.arange(9.5,12,0.1)
        bin_centers = (bins[:-1]+bins[1:])/2.0
        f_cen_red = cu.f_prop(mock['Mstar'],bins,red,blue,hosts)
        f_sat_red = cu.f_prop(mock['Mstar'],bins,red,blue,subs)
    
        p1, = plt.plot(bin_centers,f_sat_red, color='orange', alpha=alphas[i])
        p2, = plt.plot(bin_centers,f_cen_red, color='green', alpha=alphas[i])
        p1s.append(p1)
        p2s.append(p2)
    rhos = (r'$\sigma=0.0$',r'$\sigma=0.1$',r'$\sigma=0.2$',r'$\sigma=0.3$')
    plt.legend(tuple(p1s),rhos, frameon=False, fontsize=10.0, loc=4, title='satellites',labelspacing=0.01)
    plt.xlim([9.5,11.5])
    plt.ylim([0,1])
    plt.xlabel(r'$M_{\rm *} [h^{-2}M_{\odot}]$')
    plt.ylabel(r'$f_{\rm q}$')
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'sat_quenching_sigma_effect'
    fig.savefig(savepath+filename+'.pdf')
    

if __name__ == '__main__':
    main()
