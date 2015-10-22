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

    rbins = np.linspace(-1.0,1.4,25)
    rbins = 10.0**rbins
    rbins = (rbins[:-1]+rbins[1:])/2.0
    N_samples = 8*8*8

    filepath = './'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_all_9.5_10.0_result.npy'
    result_1_all = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_all_10.0_10.5_result.npy'
    result_2_all = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_all_10.5_11.0_result.npy'
    result_3_all = np.load(filepath+filename)
    
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_q_9.5_10.0_result.npy'
    result_1_q = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_q_10.0_10.5_result.npy'
    result_2_q = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_q_10.5_11.0_result.npy'
    result_3_q = np.load(filepath+filename)
    
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_sf_9.5_10.0_result.npy'
    result_1_sf = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_sf_10.0_10.5_result.npy'
    result_2_sf = np.load(filepath+filename)
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_sf_10.5_11.0_result.npy'
    result_3_sf = np.load(filepath+filename)
    
    result_all = np.hstack((result_1_all, result_2_all, result_3_all,\
                            result_1_q, result_2_q, result_3_q,\
                            result_1_sf, result_2_sf, result_3_sf))
    
    print(np.shape(result_all))
    
    all_full = result_all[0,:]
    all_sub = result_all[1:,:]
    
    cov = covariance_matrix(all_sub,all_full,N_samples)
    cov = np.matrix(cov)
    diag = np.diagonal(cov)
    
    np.save('./big_cov_8_8_8',cov)
    
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
    fig.savefig('/Users/duncan/Desktop/inv_cov_matrix.pdf', dpi=400)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.imshow(np.log10(cov), interpolation='nearest')
    cbaxes = fig.add_axes([0.87, 0.2, 0.04, 0.7]) 
    cb = plt.colorbar(cax = cbaxes)  
    plt.show()
    fig.savefig('/Users/duncan/Desktop/cov_matrix.pdf', dpi=400)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.imshow(np.fabs(corr), interpolation='nearest',vmin=0.0,vmax=1.0)
    cbaxes = fig.add_axes([0.87, 0.2, 0.04, 0.7]) 
    cb = plt.colorbar(cax = cbaxes)  
    plt.show()
    fig.savefig('/Users/duncan/Desktop/cor_matrix.pdf', dpi=400)
    
def covariance_matrix(sub,full,N_sub_vol):
    """
    Calculate the full covariance matrix.
    """
    Nr = full.shape[0] # Nr is the number of radial bins
    cov = np.zeros((Nr,Nr)) # 2D array that keeps the covariance matrix 
    after_subtraction = sub - np.mean(sub,axis=0)
    tmp = 0
    for i in range(Nr):
        for j in range(Nr):
            tmp = 0.0
            for k in range(N_sub_vol):
                tmp = tmp + after_subtraction[k,i]*after_subtraction[k,j]
            cov[i,j] = (((N_sub_vol-1)/N_sub_vol)*tmp)
    
    return cov

if __name__ == '__main__':
    main()