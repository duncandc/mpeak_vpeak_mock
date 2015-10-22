#!/usr/bin/env python

#Duncan Campbell
#February 3, 2015
#Yale University
#plot the projected correlation function of a set of mocks.

#load packages
from __future__ import print_function
import sys
import numpy as np
import h5py
from astropy.io import ascii
import matplotlib.pyplot as plt
import custom_utilities as cu

def main():

    catalogues_1 = ['sm_9.5_s0.2_sfr_c-1.0_250','sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle']

    samples = ['all','q','sf']
    sm_bins = ['9.0_9.5', '9.5_10.0', '10.0_10.5', '10.5_11.0', '11.0_11.5']

    sm_bin = sm_bins[2]

    #open correlation functions
    filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    names = ['r','wp']
    filename = catalogues_1[0]+'_wp_'+samples[0]+'_'+sm_bin+'.dat'
    result_1a = ascii.read(filepath+filename,names=names)
    filename = catalogues_1[0]+'_wp_'+samples[1]+'_'+sm_bin+'.dat'
    result_1b = ascii.read(filepath+filename,names=names)
    filename = catalogues_1[0]+'_wp_'+samples[2]+'_'+sm_bin+'.dat'
    result_1c = ascii.read(filepath+filename,names=names)
    
    filename = catalogues_1[1]+'_wp_'+samples[0]+'_'+sm_bin+'.dat'
    result_2a = ascii.read(filepath+filename,names=names)
    filename = catalogues_1[1]+'_wp_'+samples[1]+'_'+sm_bin+'.dat'
    result_2b = ascii.read(filepath+filename,names=names)
    filename = catalogues_1[1]+'_wp_'+samples[2]+'_'+sm_bin+'.dat'
    result_2c = ascii.read(filepath+filename,names=names)

    ######################################################################################
    #set up plot
    ######################################################################################
    fig1, axes = plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True,figsize=(3.3, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    
    ax = axes
    ax.set_xlim([0.1,20])
    ax.set_ylim([5,300])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\omega_p(r_p)\times r_p$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'10.0$<\log(M_{*}/M_{\odot}h^{-2})<10.5$')
    p1a, = ax.plot(result_1a['r'],result_1a['wp']*result_1a['r'],'-',color='black', alpha=1)
    p1b, = ax.plot(result_1b['r'],result_1b['wp']*result_1a['r'],'-',color='red', alpha=1)
    p1c, = ax.plot(result_1c['r'],result_1c['wp']*result_1a['r'],'-',color='blue', alpha=1)
    p2a, = ax.plot(result_2a['r'],result_2a['wp']*result_2a['r'],'--',color='black', alpha=1)
    p2b, = ax.plot(result_2b['r'],result_2b['wp']*result_2a['r'],'--',color='red', alpha=1)
    p2c, = ax.plot(result_2c['r'],result_2c['wp']*result_2a['r'],'--',color='blue', alpha=1)
    plt.show()


if __name__ == '__main__':
    main()