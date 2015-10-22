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
import sys

def main():

    sample='q'
    sm_bin='10.0_10.5'
    catalogue = 'sm_9.5_s0.2_sfr_c-0.75_250'

    #load in fiducial mock
    filepath = './'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'_cov.npy'
    cov = np.matrix(np.load(filepath+filename))
    diag = np.diagonal(cov)
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu = np.array(data['wp'])
    
    #load in comparison mock
    
    
    
    
    plt.figure()
    plt.errorbar(rbins, mu, yerr=np.sqrt(np.diagonal(cov)), color='black')
    plt.plot(rbins, wp,  color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    inv_cov = cov.I
    Y = np.matrix((wp-mu))
    
    X = Y*inv_cov*Y.T
    
    print(X)


if __name__ == '__main__':
    main()