#Duncan Campbell
#March 2015
#Yale University
#examine the relation between mpeak, vpeak, and mstar

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu

def main():

    rho = -1.0
    sigma = 0.2

    #open mock
    catalogue = 'sm_9.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_250_cen_sat_shuffle'
    #catalogue = 'sm_9.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_250_cen_shuffle'
    #catalogue = 'sm_9.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_250'
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    Lbox = 250.0
    
    hosts = (mock['upid']==-1)
    subs = (mock['upid']!=-1)
    
    #calculate the mean Vpeak Vs Mpeak
    bins = np.arange(10,15,0.2)
    bins = 10.0**bins
    bin_centers = (bins[:-1]+bins[1:])/2.0
    bins, mean_vpeak = cu.statistics.binned_mean(mock,bins=bins,\
                                                 bin_key='Mpeak',value_key='Vpeak')
    
    from scipy import interpolate
    fv_m = interpolate.InterpolatedUnivariateSpline(np.log10(bin_centers),mean_vpeak, k=1)
    
    """
    plt.figure()
    plt.plot(np.log10(bin_centers),mean_vpeak,'o')
    plt.plot(np.log10(bin_centers),fv_m(np.log10(bin_centers)))
    plt.show()
    """
    
    #subtract mean from all vpeaks
    vvpeak = mock['Vpeak'] - fv_m(np.log10(mock['Mpeak']))
    
    """
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(3.3, 3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.scatter(mock['Mpeak'],vvpeak,c=mock['SSFR'],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.', lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.ylim([-100,100])
    plt.xlim([10**10,10**15])
    plt.show()
    """
    
    #plot for a stellar mass bin
    sm_low = 10.0
    sm_high = 10.5
    selection = (mock['Mstar']>sm_low) & (mock['Mstar']<sm_high)
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(3.3, 3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.scatter(mock['Mpeak'][selection],vvpeak[selection],c=mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.', lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.ylim([-100,100])
    plt.xlim([10**10,10**15])
    plt.ylabel(r'$V_{\rm peak}-\langle V_{\rm peak} \rangle ~[{\rm km/s}]$')
    plt.xlabel(r'$M_{\rm peak} ~ [h^{-1}M_{\odot}]$')
    plt.title(str(sm_low)+r'$<\log(M_{*}/h^{-2}M_{\odot})<$'+str(sm_high))
    plt.text(10**10.1,-90,r'$\rho_{\rm SSFR}=$'+str(rho))
    plt.text(10**10.1,-75,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.show()
    
    selection = (mock['Mstar']>sm_low) & (mock['Mstar']<sm_high)
    selection = (selection & hosts)
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(3.3, 3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.scatter(mock['Mpeak'][selection],vvpeak[selection],c=mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.', lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.ylim([-100,100])
    plt.xlim([10**10,10**15])
    plt.ylabel(r'$V_{\rm peak}-\langle V_{\rm peak} \rangle ~[{\rm km/s}]$')
    plt.xlabel(r'$M_{\rm peak} ~ [h^{-1}M_{\odot}]$')
    plt.title(str(sm_low)+r'$<\log(M_{*}/h^{-2}M_{\odot})<$'+str(sm_high))
    plt.text(10**10.1,-90,r'$\rho_{\rm SSFR}=$'+str(rho))
    plt.text(10**10.1,-75,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.show(block=False)
    
    selection = (mock['Mstar']>sm_low) & (mock['Mstar']<sm_high)
    selection = (selection & subs)
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(3.3, 3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.scatter(mock['Mpeak'][selection],vvpeak[selection],c=mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.', lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.ylim([-100,100])
    plt.xlim([10**10,10**15])
    plt.ylabel(r'$V_{\rm peak}-\langle V_{\rm peak} \rangle ~[{\rm km/s}]$')
    plt.xlabel(r'$M_{\rm peak} ~ [h^{-1}M_{\odot}]$')
    plt.title(str(sm_low)+r'$<\log(M_{*}/h^{-2}M_{\odot})<$'+str(sm_high))
    plt.text(10**10.1,-90,r'$\rho_{\rm SSFR}=$'+str(rho))
    plt.text(10**10.1,-75,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.show()
    
    


if __name__ == '__main__':
    main()