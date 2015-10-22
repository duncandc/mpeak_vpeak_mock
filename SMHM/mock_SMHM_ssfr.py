#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#make plot of stellar mass vs halo mass with SSFR color-coded.

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
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle'
        #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        #catalogue = 'age_matching'
        #catalogue = 'sm_9.5_s0.2_sfr_c1.0_Halfmass_Scale_250'
        rho = -1.0
        sigma = 0.2
        
    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(catalogue)
    print(mock.dtype.names)
    
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    
    #plot stellar mass halo mass relation for the mock
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.78, bottom=0.2, top=0.9)
    p = plt.scatter(mock['mvir'][host],10.0**mock['Mstar'][host],c=mock['SSFR'][host],\
                    cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',
                    lw=0, s=2, rasterized=True)
    plt.text(10.0**10.0,5*10**11.0,r'$\rho_{\rm SSFR}=$'+str(rho))
    plt.text(10.0**10.0,5*10**10.8,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.xlabel(r'$M_{\rm vir}~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$M_{*} ~[h^{-2}M_{\odot}]$')
    plt.xlim([10.0**9.5,10.0**15.0])
    plt.ylim([10.0**9.0,10.0**12.0])
    plt.yscale('log')
    plt.xscale('log')
    plt.title("Mpeak-Vpeak mock")
    cbar_ax = fig.add_axes([0.79, 0.2, 0.02, 0.7]) #xmin, ymin, +dx, +dy
    cbar = fig.colorbar(p, cax=cbar_ax)
    cbar.set_label(r'$\log({\rm SSFR}/{\rm yr}^{-1})$')
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'SMHM'
    fig.savefig(savepath+filename+'.pdf')

if __name__ == '__main__':
    main()