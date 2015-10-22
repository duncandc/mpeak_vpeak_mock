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
from mocks.make_mock import return_mock
from halotools.mock_observables.clustering import wp

def main():

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
    
    #load covariance matrix
    filepath = './'
    filename = 'big_cov_8_8_8.npy'
    cov = np.matrix(np.load(filepath+filename))
    inv_cov = cov.I
    
    
    #mock parameters
    rhos = np.linspace(0.0,-1.0,25)
    sigmas = np.linspace(0.0,0.3,25)
    dm_sim ='Chinchilla_250'
    box_size=250.0
    period = np.array([box_size]*3)
    sm_lim = 9.5
    
    save_filename_1 = 'chinchilla_rho_3_sigma_2_jacknife_8_8_8_chi_2_array'
    save_filename_2 = 'chinchilla_rho_3_sigma_2_jacknife_8_8_8_uber_array'
    
    #correlation function parameters
    N_threads=4
    
    rp_bins = np.linspace(-1.0,1.4,25)
    rp_bins = 10.0**rp_bins
    rp_bin_centers = (rp_bins[:-1]+rp_bins[1:])/2.0
            
    pi_bins = np.linspace(0,50,50)
    pi_bin_centers = (pi_bins[:-1]+pi_bins[1:])/2.0
    
    X = np.zeros((len(sigmas),len(rhos)))
    uber_result = np.zeros((len(sigmas),len(rhos), 24*3*3))
    for i,sigma in enumerate(sigmas):
        for j,rho in enumerate(rhos):
            print(i,j,sigma, rho)
            
            mock = return_mock(sm_lim = sm_lim, sigma = sigma, rho = rho, sim = dm_sim)
            
            #calculate correlation functions
            #star forming and quenched
            LHS = -11.0
            blue  = (mock['SSFR']>LHS) #indices of blue galaxies
            red   = (mock['SSFR']<LHS) #indicies of red galaxies
            
            sm_low=9.5
            sm_high=10.0
            
            #full sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection = (selection_1 & selection_2)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of all galaxies: {0}".format(len(sample1)))
            
            result_all_a = wp(sample1, rp_bins, pi_bins, period=period,
                    do_auto=True, do_cross=False, estimator='Natural', 
                    N_threads=N_threads,max_sample_size=int(1e7))
            
            #quenched sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (red)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of red galaxies: {0}".format(len(sample1)))
            
            result_q_a = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            #star-forming sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (blue)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of blue galaxies: {0}".format(len(sample1)))
            
            result_sf_a = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            sm_low=10.0
            sm_high=10.5
            
            #full sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection = (selection_1 & selection_2)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of all galaxies: {0}".format(len(sample1)))
            
            result_all_b = wp(sample1, rp_bins, pi_bins, period=period,
                    do_auto=True, do_cross=False, estimator='Natural', 
                    N_threads=N_threads,max_sample_size=int(1e7))
            
            #quenched sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (red)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of red galaxies: {0}".format(len(sample1)))
            
            result_q_b = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            #star-forming sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (blue)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of blue galaxies: {0}".format(len(sample1)))
            
            result_sf_b = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            sm_low=10.5
            sm_high=11.0
            
            #full sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection = (selection_1 & selection_2)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of all galaxies: {0}".format(len(sample1)))
            
            result_all_c = wp(sample1, rp_bins, pi_bins, period=period,
                    do_auto=True, do_cross=False, estimator='Natural', 
                    N_threads=N_threads,max_sample_size=int(1e7))
            
            #quenched sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (red)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of red galaxies: {0}".format(len(sample1)))
            
            result_q_c = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            #star-forming sample auto correlation
            selection_1 = (mock['Mstar']>sm_low)
            selection_2 = (mock['Mstar']<sm_high)
            selection_3 = (blue)
            selection = (selection_1 & selection_2 & selection_3)
            sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
            print("number of blue galaxies: {0}".format(len(sample1)))
            
            result_sf_c = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads,max_sample_size=int(1e7))
            
            result = np.hstack((result_all_a, result_all_b, result_all_c,\
                                result_q_a, result_q_b, result_q_c,\
                                result_sf_a, result_sf_b, result_sf_c))
            
            uber_result[i,j] = result
            
            Y = np.matrix((result-mu))
            
            X_p = Y*inv_cov*Y.T
            X[i,j] = X_p
            
            print(X_p)
    
    np.save('./'+save_filename_1,X)
    np.save('./'+save_filename_2,uber_result)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.imshow(X, origin='lower',\
               extent=[sigmas[0],sigmas[-1],rhos[0],rhos[-1]],interpolation='nearest',\
               aspect='auto',cmap='cool')
    plt.plot([0.2],[-0.8],'x',color='black',ms=10)
    plt.xlabel(r'$\sigma_{\rm SMHM}$')
    plt.ylabel(r'$\rho_{\rm SSFR}$')
    plt.colorbar()
    plt.show(block=False)
    fig.savefig('/Users/duncan/Desktop/chi_squared_plot.pdf')

if __name__ == '__main__':
    main()