#!/usr/bin/env python

#Duncan Campbell
#February 3, 2015
#Yale University
#plot delta sigma of a set of mocks and compare to the SDSS values

#load packages
from __future__ import division, print_function
import sys
import numpy as np
import h5py
from astropy.io import ascii
import matplotlib.pyplot as plt
import custom_utilities as cu

def main():

    #open SDSS results
    #read in SDSS projected two-point correlation function
    names = ['r','sm9.8','sm9.8err','sm10.2','sm10.2err','sm10.6','sm10.6err']
    sdss_delta_sigma_sf= ascii.read("./Watson_2015_delta_sigma_sf.txt",names=names) 
    sdss_delta_sigma_q= ascii.read("./Watson_2015_delta_sigma_q.txt",names=names)
    sdss_delta_sigma_all= ascii.read("./Hearin_2014_delta_sigma_all.txt",names=names)
    h=0.7 #for the SDSS measurements, make sure to scale correctly when plotting
    
    catalogues = ['sm_9.49_s0.2_sfr_c-1.0_250','sm_9.49_s0.2_sfr_c-0.75_250',\
                  'sm_9.49_s0.2_sfr_c-0.5_250','sm_9.49_s0.2_sfr_c-0.25_250',\
                  'sm_9.49_s0.2_sfr_c0.0_250']
    samples = ['all','q','sf']
    sm_thresholds = ['9.49','9.89','10.29']
    
    fig1, axes = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True,figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05, left=0.1, right=0.95, bottom=0.2, top=0.9)
    
    r_all = sdss_delta_sigma_all['r']*h/1000.0
    r_q = sdss_delta_sigma_q['r']*h/1000.0
    r_sf = sdss_delta_sigma_sf['r']*h/1000.0
    
    ax = axes[0]
    #ax.errorbar(r_all,sdss_delta_sigma_all['sm9.8']/h,\
    #            yerr=sdss_delta_sigma_all['sm9.8err']/h,\
    #            ms=3, fmt='o', color='black', mec='none')
    #ax.errorbar(r_q,sdss_delta_sigma_q['sm9.8']/h,\
    #            yerr=sdss_delta_sigma_q['sm9.8err']/h,\
    #            ms=3, fmt='o',color='red', mec='none')
    #ax.errorbar(r_sf,sdss_delta_sigma_sf['sm9.8']/h,\
    #            yerr=sdss_delta_sigma_sf['sm9.8err']/h,\
    #            ms=3, fmt='o',color='blue', mec='none')
    ax.set_xlim([0.1,2])
    ax.set_ylim([1,100])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\Delta\Sigma(r_p)~[M_{\odot}{\rm pc}^{-2}h]$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$\log(M_{*}~M_{\odot}h^{-2})>9.49$')
    
    ax = axes[1]
    #ax.errorbar(r_all,sdss_delta_sigma_all['sm10.2']/h,\
    #            yerr=sdss_delta_sigma_all['sm10.2err']/h,\
    #            ms=3, fmt='o', color='black', mec='none')
    #ax.errorbar(r_q,sdss_delta_sigma_q['sm10.2']/h,\
    #            yerr=sdss_delta_sigma_q['sm10.2err']/h,\
    #            ms=3, fmt='o',color='red', mec='none')
    #ax.errorbar(r_sf,sdss_delta_sigma_sf['sm10.2']/h,\
    #            yerr=sdss_delta_sigma_sf['sm10.2err']/h,\
    #            ms=3, fmt='o',color='blue', mec='none')
    ax.set_xlim([0.1,2])
    ax.set_ylim([1,100])
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylabel(r'$\omega_p(r_p)$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$\log(M_{*}~M_{\odot}h^{-2})>9.89$')
    
    ax = axes[2]
    #ax.errorbar(r_all,sdss_delta_sigma_all['sm10.6']/h,\
    #            yerr=sdss_delta_sigma_all['sm10.6err']/h,\
    #            ms=3, fmt='o', color='black', mec='none')
    #ax.errorbar(r_q,sdss_delta_sigma_q['sm10.6']/h,\
    #            yerr=sdss_delta_sigma_q['sm10.6err']/h,\
    #            ms=3, fmt='o',color='red', mec='none')
    #ax.errorbar(r_sf,sdss_delta_sigma_sf['sm10.6']/h,\
    #            yerr=sdss_delta_sigma_sf['sm10.6err']/h,\
    #            ms=3, fmt='o',color='blue', mec='none')
    ax.set_xlim([0.1,2])
    ax.set_ylim([1,100])
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_ylabel(r'$\omega_p(r_p)$')
    ax.set_xlabel(r'$r_p~[{\rm Mpc}~h^{-1}]$')
    ax.set_title(r'$\log(M_{*}~M_{\odot}h^{-2})>10.29$')
    
    for i in range(0,len(sm_thresholds)):
        filepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
        names = ['r','delta_sigma']
        filename = catalogues[0]+'_DeltaSigma_'+samples[0]+'_'+sm_thresholds[i]+'.dat'
        result_1a = ascii.read(filepath+filename,names=names)
        filename = catalogues[0]+'_DeltaSigma_'+samples[1]+'_'+sm_thresholds[i]+'.dat'
        result_1b = ascii.read(filepath+filename,names=names)
        filename = catalogues[0]+'_DeltaSigma_'+samples[2]+'_'+sm_thresholds[i]+'.dat'
        result_1c = ascii.read(filepath+filename,names=names)
    
        filename = catalogues[1]+'_DeltaSigma_'+samples[0]+'_'+sm_thresholds[i]+'.dat'
        result_2a = ascii.read(filepath+filename,names=names)
        filename = catalogues[1]+'_DeltaSigma_'+samples[1]+'_'+sm_thresholds[i]+'.dat'
        result_2b = ascii.read(filepath+filename,names=names)
        filename = catalogues[1]+'_DeltaSigma_'+samples[2]+'_'+sm_thresholds[i]+'.dat'
        result_2c = ascii.read(filepath+filename,names=names)
    
        filename = catalogues[2]+'_DeltaSigma_'+samples[0]+'_'+sm_thresholds[i]+'.dat'
        result_3a = ascii.read(filepath+filename,names=names)
        filename = catalogues[2]+'_DeltaSigma_'+samples[1]+'_'+sm_thresholds[i]+'.dat'
        result_3b = ascii.read(filepath+filename,names=names)
        filename = catalogues[2]+'_DeltaSigma_'+samples[2]+'_'+sm_thresholds[i]+'.dat'
        result_3c = ascii.read(filepath+filename,names=names)
        
        filename = catalogues[3]+'_DeltaSigma_'+samples[0]+'_'+sm_thresholds[i]+'.dat'
        result_4a = ascii.read(filepath+filename,names=names)
        filename = catalogues[3]+'_DeltaSigma_'+samples[1]+'_'+sm_thresholds[i]+'.dat'
        result_4b = ascii.read(filepath+filename,names=names)
        filename = catalogues[3]+'_DeltaSigma_'+samples[2]+'_'+sm_thresholds[i]+'.dat'
        result_4c = ascii.read(filepath+filename,names=names)
        
        filename = catalogues[4]+'_DeltaSigma_'+samples[0]+'_'+sm_thresholds[i]+'.dat'
        result_5a = ascii.read(filepath+filename,names=names)
        filename = catalogues[4]+'_DeltaSigma_'+samples[1]+'_'+sm_thresholds[i]+'.dat'
        result_5b = ascii.read(filepath+filename,names=names)
        filename = catalogues[4]+'_DeltaSigma_'+samples[2]+'_'+sm_thresholds[i]+'.dat'
        result_5c = ascii.read(filepath+filename,names=names)
    
        ax = axes[i]
        p1a, = ax.plot(result_1a['r'],result_1a['delta_sigma'],'-',color='black', alpha=1)
        p1b, = ax.plot(result_1b['r'],result_1b['delta_sigma'],'-',color='red', alpha=1)
        p1c, = ax.plot(result_1c['r'],result_1c['delta_sigma'],'-',color='blue', alpha=1)
    
        p2a, = ax.plot(result_2a['r'],result_2a['delta_sigma'],'-',color='black', alpha=0.8)
        p2b, = ax.plot(result_2b['r'],result_2b['delta_sigma'],'-',color='red', alpha=0.8)
        p2c, = ax.plot(result_2c['r'],result_2c['delta_sigma'],'-',color='blue', alpha=0.8)
    
        p3a, = ax.plot(result_3a['r'],result_3a['delta_sigma'],'-',color='black', alpha=0.6)
        p3b, = ax.plot(result_3b['r'],result_3b['delta_sigma'],'-',color='red', alpha=0.6)
        p3c, = ax.plot(result_3c['r'],result_3c['delta_sigma'],'-',color='blue', alpha=0.6)
        
        p4a, = ax.plot(result_4a['r'],result_4a['delta_sigma'],'-',color='black', alpha=0.4)
        p4b, = ax.plot(result_4b['r'],result_4b['delta_sigma'],'-',color='red', alpha=0.4)
        p4c, = ax.plot(result_4c['r'],result_4c['delta_sigma'],'-',color='blue', alpha=0.4)
        
        p5a, = ax.plot(result_5a['r'],result_5a['delta_sigma'],'-',color='black', alpha=0.2)
        p5b, = ax.plot(result_5b['r'],result_5b['delta_sigma'],'-',color='red', alpha=0.2)
        p5c, = ax.plot(result_5c['r'],result_5c['delta_sigma'],'-',color='blue', alpha=0.2)
        
        if i==0:
            ax.legend((p1a,p2a,p3a,p4a,p5a),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=1, fontsize=10, frameon=False, labelspacing=0.01, title='all')
        if i==1:
            ax.legend((p1b,p2b,p3b,p4b,p5b),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=1, fontsize=10, frameon=False, labelspacing=0.01, title='quenched')
        if i==2:
            ax.legend((p1c,p2c,p3c,p4c,p5c),\
            (r"$\rho=-1.0$",r"$\rho=-0.75$",r"$\rho=-0.5$",r"$\rho=-0.25$",r"$\rho=0.0$"),\
            loc=1, fontsize=10, frameon=False, labelspacing=0.01, title='star-forming')
    
    plt.show()
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'delta_sigma'
    fig1.savefig(savepath+filename+'.pdf')
    
    
if __name__ == '__main__':
    main()