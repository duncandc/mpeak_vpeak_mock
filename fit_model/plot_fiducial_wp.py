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

    N_top = 20 #top number of fits to show

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
    
    rp_bins = np.linspace(-1.0,1.4,25)
    rp_bins = 10.0**rp_bins
    r = (rp_bins[:-1]+rp_bins[1:])/2.0
    
    #load covariance matrix
    filepath = './'
    filename = 'big_cov_8_8_8.npy'
    cov = np.matrix(np.load(filepath+filename))
    diag = np.sqrt(np.diag(cov))
    inds = np.arange(0,24*3*3)
    a,b,c,d,e,f,g,h,i = np.split(inds,9)
    
    all_a_err = diag[a]
    all_b_err = diag[b]
    all_c_err = diag[c]
    
    red_a_err = diag[d]
    red_b_err = diag[e]
    red_c_err = diag[f]
    
    blue_a_err = diag[g]
    blue_b_err = diag[h]
    blue_c_err = diag[i]
    
    
    fig1, axes = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True,figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05, left=0.1, right=0.95, bottom=0.2, top=0.9)
    
    ax = axes[0]
    ax.set_xlim([0.1,20])
    ax.set_ylim([5,300])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$r_p\times\omega_p(r_p)$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$9.5<\log(M_{*}/M_{\odot}h^{-2})<10.0$')
    
    ax = axes[1]
    ax.set_xlim([0.1,20])
    ax.set_ylim([5,300])
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylabel(r'$\omega_p(r_p)$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$10.0<\log(M_{*}/M_{\odot}h^{-2})<10.5$')
    
    ax = axes[2]
    ax.set_xlim([0.1,20])
    ax.set_ylim([5,300])
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylabel(r'$\omega_p(r_p)$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$10.5<\log(M_{*}/M_{\odot}h^{-2})<11.0$')
    
    ax = axes[0]
    p1a = ax.errorbar(r,mu_all_a*r, yerr=all_a_err*r, fmt = 'o', color='black', ms=4, mec='none')
    p1b = ax.errorbar(r,mu_q_a*r, yerr=red_a_err*r, fmt = 'o',color='red', ms=4, mec='none')
    p1c = ax.errorbar(r,mu_sf_a*r, yerr=blue_a_err*r, fmt = 'o',color='blue', ms=4, mec='none')
    
    ax = axes[1]
    p2a = ax.errorbar(r,mu_all_b*r, yerr=all_b_err*r, fmt = 'o',color='black', ms=4, mec='none')
    p2b = ax.errorbar(r,mu_q_b*r, yerr=red_b_err*r, fmt = 'o',color='red', ms=4, mec='none')
    p2c = ax.errorbar(r,mu_sf_b*r, yerr=blue_b_err*r, fmt = 'o',color='blue', ms=4, mec='none')
    
    ax = axes[2]
    p2a = ax.errorbar(r,mu_all_c*r, yerr=all_c_err*r, fmt = 'o', color='black', ms=4, mec='none')
    p2b = ax.errorbar(r,mu_q_c*r, yerr=red_c_err*r, fmt = 'o', color='red', ms=4, mec='none')
    p2c = ax.errorbar(r,mu_sf_c*r, yerr=blue_c_err*r, fmt = 'o', color='blue', ms=4, mec='none')
    
    #open fitted models and overplot
    #load model correlation functions
    filepath = './'
    filename = 'chinchilla_rho_3_sigma_2_jacknife_8_8_8_uber_array.npy'
    uber_result = np.load(filepath+filename)
    filename = 'chinchilla_rho_3_sigma_2_jacknife_8_8_8_chi_2_array.npy'
    X = np.load(filepath+filename)
    
    inds = np.argsort(X, axis=None)
    inds = np.unravel_index(inds,(25,25))
    
    rhos = np.linspace(0.0,-1.0,25)
    sigmas = np.linspace(0.0,0.3,25)
    xx, yy = np.meshgrid(sigmas, rhos)
    
    ii = 19
    jj = 16
    print(rhos[ii],sigmas[jj])
    
    for i in range(0,N_top):
        
        print(i)
        a,b,c,d,e,f,g,h,i = np.split(uber_result[inds[0][i],inds[1][i]], 9)
    
        ax = axes[0]
        p1aa, = ax.plot(r, a*r, color='black', alpha=0.2)
        p1bb, = ax.plot(r, d*r, color='red', alpha=0.2)
        p1cc, = ax.plot(r, g*r, color='blue', alpha=0.2)
        
        ax = axes[1]
        p2aa, = ax.plot(r, b*r, color='black', alpha=0.2)
        p2bb, = ax.plot(r, e*r, color='red', alpha=0.2)
        p2cc, = ax.plot(r, h*r, color='blue', alpha=0.2)
    
        ax = axes[2]
        p2aa, = ax.plot(r, c*r, color='black', alpha=0.2)
        p2bb, = ax.plot(r, f*r, color='red', alpha=0.2)
        p2cc, = ax.plot(r, i*r, color='blue', alpha=0.2)
    
    a,b,c,d,e,f,g,h,i = np.split(uber_result[jj,ii], 9) # input model
    
    ax = axes[0]
    p1aaa, = ax.plot(r, a*r, color='black', alpha=1)
    p1bbb, = ax.plot(r, d*r, color='red', alpha=1)
    p1ccc, = ax.plot(r, g*r, color='blue', alpha=1)
        
    ax = axes[1]
    p2aaa, = ax.plot(r, b*r, color='black', alpha=1)
    p2bbb, = ax.plot(r, e*r, color='red', alpha=1)
    p2ccc, = ax.plot(r, h*r, color='blue', alpha=1)
    
    ax = axes[2]
    p2aaa, = ax.plot(r, c*r, color='black', alpha=1)
    p2bbb, = ax.plot(r, f*r, color='red', alpha=1)
    p2ccc, = ax.plot(r, i*r, color='blue', alpha=1)
    
    ax = axes[0]
    ax.legend([p1a,p1b,p1c], ('all','quenched','star-forming'), title='fiducial mock',\
              numpoints=1, fontsize=10, frameon=False, loc=4, labelspacing=0.3)
    
    ax = axes[1]
    ax.legend([p1aaa,p1aa], ['correct model', 'top 20 models'],\
              fontsize=10, frameon=False, loc=4, labelspacing=0.2)
    
    
    plt.show()
    fig1.savefig('/Users/duncan/Desktop/wp_models.pdf', dpi=300)
    
    colors = X[inds[0][0:N_top],inds[1][0:N_top]]/(24*3*3)
    
    fig2 = plt.figure(figsize=(3.3,3.3))
    fig2.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.scatter(xx.T[inds[0][0:N_top],inds[1][0:N_top]], yy.T[inds[0][0:N_top],inds[1][0:N_top]],\
             c = colors, lw=0, cmap='cool', marker='.')
    plt.xlim([0,0.3])
    plt.ylim([0,-1])
    plt.show(block=False)
    
if __name__ == '__main__':
    main()