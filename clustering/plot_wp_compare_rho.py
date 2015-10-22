#!/usr/bin/env python

#Duncan Campbell
#February 3, 2015
#Yale University
#plot the projected correlation function of a set of mocks and compare to the SDSS values

#load packages
from __future__ import print_function
import sys
import numpy as np
import h5py
from astropy.io import ascii
import matplotlib.pyplot as plt
import custom_utilities as cu

def main():

    catalogues_0 = ['sm_8.5_s0.2_sfr_c-1.0_250','sm_8.5_s0.2_sfr_c-0.75_250',\
                  'sm_8.5_s0.2_sfr_c-0.5_250','sm_8.5_s0.2_sfr_c-0.25_250',\
                  'sm_8.5_s0.2_sfr_c0.0_250']
    catalogues_1 = ['sm_9.5_s0.2_sfr_c-1.0_250','sm_9.5_s0.2_sfr_c-0.75_250',\
                  'sm_9.5_s0.2_sfr_c-0.5_250','sm_9.5_s0.2_sfr_c-0.25_250',\
                  'sm_9.5_s0.2_sfr_c0.0_250']
    samples = ['all','q','sf']
    sm_bins = ['9.0_9.5', '9.5_10.0', '10.0_10.5', '10.5_11.0', '11.0_11.5']

    ######################################################################################
    #set up plot
    ######################################################################################
    fig1, axes = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True,figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05, left=0.1, right=0.95, bottom=0.2, top=0.9)
    
    ax = axes[0]
    ax.set_xlim([0.1,20])
    ax.set_ylim([5,300])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\omega_p(r_p)\times r_p$')
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

    ######################################################################################
    #read in results and plot
    ######################################################################################

    for i,sm_bin in enumerate(sm_bins[1:4]):
        
        print(i, sm_bin)
        
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
    
        filename = catalogues_1[2]+'_wp_'+samples[0]+'_'+sm_bin+'.dat'
        result_3a = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[2]+'_wp_'+samples[1]+'_'+sm_bin+'.dat'
        result_3b = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[2]+'_wp_'+samples[2]+'_'+sm_bin+'.dat'
        result_3c = ascii.read(filepath+filename,names=names)
    
        filename = catalogues_1[3]+'_wp_'+samples[0]+'_'+sm_bin+'.dat'
        result_4a = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[3]+'_wp_'+samples[1]+'_'+sm_bin+'.dat'
        result_4b = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[3]+'_wp_'+samples[2]+'_'+sm_bin+'.dat'
        result_4c = ascii.read(filepath+filename,names=names)
    
        filename = catalogues_1[4]+'_wp_'+samples[0]+'_'+sm_bin+'.dat'
        result_5a = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[4]+'_wp_'+samples[1]+'_'+sm_bin+'.dat'
        result_5b = ascii.read(filepath+filename,names=names)
        filename = catalogues_1[4]+'_wp_'+samples[2]+'_'+sm_bin+'.dat'
        result_5c = ascii.read(filepath+filename,names=names)
    
        ax = axes[i]
        
        p1a, = ax.plot(result_1a['r'],result_1a['wp']*result_1a['r'],'-',color='black', alpha=1)
        p1b, = ax.plot(result_1b['r'],result_1b['wp']*result_1a['r'],'-',color='red', alpha=1)
        p1c, = ax.plot(result_1c['r'],result_1c['wp']*result_1a['r'],'-',color='blue', alpha=1)
    
        p2a, = ax.plot(result_2a['r'],result_2a['wp']*result_1a['r'],'-',color='black', alpha=0.8)
        p2b, = ax.plot(result_2b['r'],result_2b['wp']*result_1a['r'],'-',color='red',alpha=0.8)
        p2c, = ax.plot(result_2c['r'],result_2c['wp']*result_1a['r'],'-',color='blue',alpha=0.8)
    
        p3a, = ax.plot(result_3a['r'],result_3a['wp']*result_1a['r'],'-',color='black',alpha=0.6)
        p3b, = ax.plot(result_3b['r'],result_3b['wp']*result_1a['r'],'-',color='red', alpha=0.6)
        p3c, = ax.plot(result_3c['r'],result_3c['wp']*result_1a['r'],'-',color='blue', alpha=0.6)
    
        p4a, = ax.plot(result_4a['r'],result_4a['wp']*result_1a['r'],'-',color='black',alpha=0.4)
        p4b, = ax.plot(result_4b['r'],result_4b['wp']*result_1a['r'],'-',color='red', alpha=0.4)
        p4c, = ax.plot(result_4c['r'],result_4c['wp']*result_1a['r'],'-',color='blue', alpha=0.4)
    
        p5a, = ax.plot(result_5a['r'],result_5a['wp']*result_1a['r'],'-',color='black',alpha=0.2)
        p5b, = ax.plot(result_5b['r'],result_5b['wp']*result_1a['r'],'-',color='red', alpha=0.2)
        p5c, = ax.plot(result_5c['r'],result_5c['wp']*result_1a['r'],'-',color='blue', alpha=0.2)
    
        if i==0:
            ax.legend((p1a,p2a,p3a,p4a,p5a),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=4,fontsize=10, frameon=False, labelspacing=0.01, title='all')
        if i==1:
            ax.legend((p1b,p2b,p3b,p4b,p5b),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=4,fontsize=10, frameon=False, labelspacing=0.01, title='quenched')
        if i==2:
            ax.legend((p1c,p2c,p3c,p4c,p5c),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=4,fontsize=10, frameon=False, labelspacing=0.01, title='star-forming')
        
    plt.show(block=True)
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'w_p_cross_r'
    fig1.savefig(savepath+filename+'.pdf')


if __name__ == '__main__':
    main()