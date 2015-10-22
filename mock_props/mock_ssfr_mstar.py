#!/usr/bin/env python

#Author: Duncan Campbell
#January 28, 2015
#Yale University
#plot the SSFR vs stellar for mock galaxies

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from astropy.cosmology import FlatLambdaCDM

def main():
    
    
    catalogue_1 = 'sm_9.5_s0.2_sfr_c0.0_250'
    catalogue_2 = 'sm_8.5_s0.2_sfr_c-1.0_125'
    
    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
    #open catalogue
    f_1 = h5py.File(filepath_mock+catalogue_1+'.hdf5', 'r') #open catalogue file
    mock_1 = f_1.get(catalogue_1)
    f_2 = h5py.File(filepath_mock+catalogue_2+'.hdf5', 'r') #open catalogue file
    mock_2 = f_2.get(catalogue_2)
    data_1 = np.array(mock_1)
    data_2 = np.array(mock_2)
    
    GC = get_sdss_sample()
    
    mock = np.concatenate((data_1,data_2),axis=1)
    
    sm_bins = np.arange(9.5,11.5,0.4)
    sm_bin_centers = (sm_bins[:-1]+sm_bins[1:])/2.0
    binned_gc_inds = np.digitize(mock['Mstar'],bins=sm_bins)
    p = []
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    ax = fig.add_subplot(1, 1, 1)
    for i in range(0,len(sm_bins)-1):
        gc_inds = np.where(binned_gc_inds==i+1)[0]
        sub_sfr = mock['SSFR'][gc_inds]
        N = len(gc_inds)
        sorted_sfr = np.sort(sub_sfr)
        count = np.ones(N)
        counts = np.cumsum(count)-1
        normalized_counts = counts/float(counts[-1])
        
        line, = ax.plot(sorted_sfr,normalized_counts,color='black',alpha=(i+1.0)/(len(sm_bins)))
        p.append(line)
    
    
    binned_gc_inds = np.digitize(np.log10(10.0**GC['sm_MEDIAN']),bins=sm_bins)
    for i in range(0,len(sm_bins)-1):
        gc_inds = np.where(binned_gc_inds==i+1)[0]
        sub_sfr = GC['sfr_MEDIAN'][gc_inds]
        N = len(gc_inds)
        sorted_sfr = np.sort(sub_sfr)
        count = np.ones(N)
        counts = np.cumsum(count)-1
        normalized_counts = counts/float(counts[-1])
        
        ax.plot(sorted_sfr,normalized_counts,'--',color='black',alpha=(i+1.0)/(len(sm_bins)))
    
    ax.legend(p,sm_bin_centers,loc=4,title=r'$M_{*}$',fontsize=10, frameon=False, labelspacing=0.01)
    ax.set_ylim([0,1])
    ax.set_xlim([-12.5,-8.5])
    ax.set_xlabel(r'$\log({\rm SSFR}/{\rm yr^{-1}})$')
    ax.set_ylabel(r'$f(<{\rm SSFR})$')
    ax.set_xticklabels(['-12.5','','-11.5','','-10.5','','-9.5','','-8.5'])
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'mstar_ssfr_dist'
    fig.savefig(savepath+filename+'.pdf')


def get_sdss_sample():
    
    from scipy import interpolate
    
    filepath = cu.get_output_path() + 'processed_data/NYU_VAGC/'
    galaxy_catalogue = 'nyu_lss_mpa_vagc_dr7'
    f =  h5py.File(filepath+galaxy_catalogue+'.hdf5', 'r')
    GC = f.get(galaxy_catalogue) #halo catalogue
    GC = np.array(GC)
    for name in GC.dtype.names: print(name)
    
    #trim the catalogue a bit
    zmin = 0.01
    zmax = 0.2
    selection = (GC['M']<=17.6) & (GC['Z']<zmax) & (GC['Z']>zmin) &\
                (GC['ZTYPE']==1) & (GC['FGOTMAIN']>0)
    GC = GC[selection]
    
    sm_key = 'sm_MEDIAN'
    GC[sm_key] = GC[sm_key]+np.log10(0.7**2.0)
    
    #make cuts on data to get a clean and complete sample
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    z = np.linspace(0.001,1,1000)
    dL = cosmo.luminosity_distance(z).value
    #make cheater fit to dl(z) function
    dL_z = interpolate.interp1d(z,dL,kind='linear')
    Mstar_lim = (4.852 + 2.246*np.log10(dL) + 1.123*np.log10(1+z) - 1.186*z)/(1.0-0.067*z)
    #make completeness cut
    #dL = cosmo.luminosity_distance(GC['Z']).value
    dL = dL_z(GC['Z'])
    z = GC['Z']
    LHS = (4.852 + 2.246*np.log10(dL) + 1.123*np.log10(1.0+z) - 1.186*z)/(1.0-0.067*z)
    keep = (GC[sm_key]>LHS)
    GC = GC[keep]
    
    return GC

if __name__ == '__main__':
    main()