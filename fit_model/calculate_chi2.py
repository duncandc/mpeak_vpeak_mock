#!/usr/bin/env python

#Duncan Campbell
#September 2015
#Yale University
#from a grid of parameters, make mocks

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.io import ascii
import custom_utilities as cu
import sys
from scipy.ndimage.filters import gaussian_filter

def main():

    #load covariance matrix
    filepath = './'
    filename = 'big_cov_8_8_8.npy'
    cov = np.matrix(np.load(filepath+filename))
    inv_cov = cov.I
    
    cov_p = np.matrix(np.zeros((np.shape(cov)[0],np.shape(cov)[1])))
    np.fill_diagonal(cov_p,np.diagonal(cov))
    
    cov = cov_p
    inv_cov = cov.I
    
    #load in fiducial mock
    #all galaxies
    sample='all'
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    sm_bin = '9.5_10.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_all_a = np.array(data['wp'])
    sm_bin = '10.0_10.5'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_all_b = np.array(data['wp'])
    sm_bin = '10.5_11.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_all_c = np.array(data['wp'])
    #quenched galaxies
    sample='q'
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    sm_bin = '9.5_10.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_q_a = np.array(data['wp'])
    sm_bin = '10.0_10.5'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_q_b = np.array(data['wp'])
    sm_bin = '10.5_11.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_q_c = np.array(data['wp'])
    #SF galaxies
    sample='sf'
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    sm_bin = '9.5_10.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_sf_a = np.array(data['wp'])
    sm_bin = '10.0_10.5'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_sf_b = np.array(data['wp'])
    sm_bin = '10.5_11.0'
    filename = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250_wp_fiducial_'+sample+'_'+sm_bin+'.dat'
    data = ascii.read(filepath+filename)
    rbins = np.array(data['r'])
    mu_sf_c = np.array(data['wp'])

    mu = np.hstack((mu_all_a, mu_all_b, mu_all_c,mu_q_a, mu_q_b, mu_q_c,mu_sf_a, mu_sf_b, mu_sf_c))
    
    #load model correlation functions
    filepath = './'
    filename = 'rho_25_sigma_25_jacknife_8_8_8_uber_array.npy'
    uber_result = np.load(filepath+filename)
    
    #parameters
    rhos = np.linspace(0.0,-1.0,25)
    sigmas = np.linspace(0.0,0.3,25)
    
    X = np.zeros((len(sigmas),len(rhos)))
    for i,sigma in enumerate(sigmas):
        for j,rho in enumerate(rhos):
            print(i,j,sigma, rho)
            
            Y = np.matrix((uber_result[i,j]-mu))
            
            X_p = Y*inv_cov*Y.T
            X[i,j] = X_p
    
    X = X/(24*3*3)
    min_mask = (X==np.min(X))
    
    print(np.min(X))
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.imshow(X.T, origin='lower',\
               extent=[sigmas[0],sigmas[-1],rhos[0],rhos[-1]],interpolation='nearest',\
               aspect='auto',cmap='hot', vmin=0.0)
    cbar = plt.colorbar()
    cbar.set_label(r'$\chi^2$')
    xx, yy = np.meshgrid(sigmas, rhos)
    CS = plt.contour(xx, yy, gaussian_filter(X.T,(1,1)), colors='white')
    plt.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
    plt.plot([0.2],[-0.8],'x',color='white',ms=10)
    plt.plot([xx[min_mask.T]],[yy[min_mask.T]],'x',color='grey',ms=10)
    plt.xlabel(r'$\sigma_{\rm SMHM}$')
    plt.ylabel(r'$\rho_{\rm SSFR}$')
    plt.show()
    
if __name__ == "__main__":
    main()