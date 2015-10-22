#!/usr/bin/env python

#Duncan Campbell
#September 2015
#Yale University
#examine the covariance matrix of the correlation function for our fiducial mock

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.io import ascii
import custom_utilities as cu

def main():

    sample='9.5_10.0'

    filepath = './'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_all_'+sample+'_cov.npy'
    cov = np.load(filepath+filename)
    cov = np.matrix(cov)
    diag = np.diagonal(cov)
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_all_'+sample+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    rbins = np.log10(rbins)
    mu = np.array(data['wp'])
    
    corr = np.zeros(np.shape(cov))
    for i in range(0,np.shape(cov)[0]):
        for j in range(0,np.shape(cov)[1]):
            corr[i,j] =  cov[i,j]/(np.sqrt(diag[i]*diag[j]))
    
    i_cov = cov.I
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.imshow(np.log10(np.fabs(i_cov)), interpolation='nearest')
    cbaxes = fig.add_axes([0.87, 0.2, 0.04, 0.7]) 
    cb = plt.colorbar(cax = cbaxes)  
    plt.show()
    #fig.savefig('/Users/duncan/Desktop/inv_cov_matrix.pdf', dpi=400)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.imshow(np.fabs(corr), interpolation='nearest',vmin=0.0,vmax=1.0,\
               extent=[rbins[0],rbins[-1],rbins[-1],rbins[0]])
    plt.xlabel(r'$\log(r_p)$')
    plt.ylabel(r'$\log(r_p)$')
    cbaxes = fig.add_axes([0.87, 0.2, 0.04, 0.7]) 
    cb = plt.colorbar(cax = cbaxes)  
    plt.show()
    #fig.savefig('/Users/duncan/Desktop/cor_matrix.pdf')



if __name__ == '__main__':
    main()