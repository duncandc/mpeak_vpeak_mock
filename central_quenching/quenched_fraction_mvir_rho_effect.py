#!/usr/bin/env python

#Author: Duncan Campbell
#January 28, 2015
#Yale University
#plot the central quenched fraction

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
from custom_utilities.cython_utilities.match2 import match2
import sys

def main():
    
    rhos = np.array([-1.0,-0.75,-0.5,-0.25,0.0])
    sigma = 0.2
    
    #fig = plt.figure(figsize=(3.3,3.3))
    #fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    #ax = fig.add_subplot(1,1,1)
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(3.3,5.5))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.1, top=0.95)
    axes = axes.flatten()
    colors = [0.8,0.6,0.4,0.2,0.0]
    colors = colors[::-1]
    pa = []
    pb = []
    for i,rho in enumerate(rhos):
        catalogue_1 = 'sm_9.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_250'
        catalogue_2 = 'sm_8.5_s'+str(sigma)+'_sfr_c'+str(rho)+'_125'
    
        filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
        print('opening mock catalogue:', catalogue_1+'.hdf5')
        #open catalogue
        f = h5py.File(filepath_mock+catalogue_1+'.hdf5', 'r') #open catalogue file
        mock_1 = f.get(catalogue_1)
        print(mock_1.dtype.names)
        print('opening mock catalogue:', catalogue_2+'.hdf5')
        #open catalogue
        f = h5py.File(filepath_mock+catalogue_2+'.hdf5', 'r') #open catalogue file
        mock_2 = f.get(catalogue_2)
    
        am_1 = abundance_match(mock_1)
        am_2 = abundance_match(mock_2)
        am = np.hstack((am_1,am_2))
    
        mock = np.hstack((mock_1,mock_2))
    
        bins = np.arange(10,15,0.1)
        bin_centers = (bins[:-1]+bins[1:])/2.0
        mass = np.log10(mock['mvir'])
        f_red_cen_halo, N_halo = quenched_fraction(mock,mass,bins)
        mass = np.log10(am)
        f_red_cen_am, N_am = quenched_fraction(mock,mass,bins,use_rank=False)
    
        #pa_1, = axes[0].plot(10**bin_centers,f_red_cen_halo,'-',color=str(colors[i]))
        pa_1 = cu.plotting.colorline(10**bin_centers, f_red_cen_halo, rho, axes[0],\
                                     vmin=np.min(rhos), vmax=np.max(rhos), cmap='winter')
        #pb_1, = axes[1].plot(10**bin_centers,f_red_cen_am,'-',color=str(colors[i]))
        pb_1 = cu.plotting.colorline(10**bin_centers, f_red_cen_am, rho, axes[1],\
                                     vmin=np.min(rhos), vmax=np.max(rhos), cmap='winter')
        pa.append(pa_1)
        pb.append(pb_1)

    axes[0].set_xlabel(r'$M_{\rm vir}$ $[h^{-1}M_{\odot}]$')
    axes[0].set_ylabel(r'$f_{\rm q}(M)$')
    axes[0].set_xlim([10.0**11,10.0**14.0])
    axes[0].set_xscale('log')
    #axes[0].set_xticklabels([' ',r'$10^{11}$',r'$10^{12}$',r'$10^{13}$',' '])
    axes[0].legend(pa,rhos,loc=4,fontsize=10,frameon=False,labelspacing=0.01,title=r'$\rho_{\rm SSFR}$')
    axes[0].text(10.0**11.2,0.9,"intrinisic")
    axes[0].text(10.0**11.2,0.84,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    
    axes[1].set_xlabel(r'$M_{\rm vir}$ $[h^{-1}M_{\odot}]$')
    axes[1].set_ylabel(r'$f_{\rm q}(M)$')
    axes[1].set_xlim([10.0**11,10.0**14.0])
    axes[1].set_xscale('log')
    #axes[1].set_xticklabels([' ',r'$10^{11}$',r'$10^{12}$',r'$10^{13}$','$10^{14}$'])
    axes[1].set_yticklabels(['0',r'$0.2$',r'$0.4$',r'$0.6$','$0.8$',' '])
    axes[1].legend(pb,rhos,loc=4,fontsize=10,frameon=False,labelspacing=0.01,title=r'$\rho_{\rm SSFR}$')
    axes[1].text(10.0**11.2,0.85,"inverse abundance \n matching")
    axes[1].text(10.0**11.2,0.78,r'$\sigma_{\rm SMHM}=$'+str(sigma))
    plt.show()
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'f_q_cen_compare_rho'
    fig.savefig(savepath+filename+'.pdf')


def quenched_fraction(mock,mass,bins,use_rank=False):
    
    if use_rank==True:
        #rank by stellar mass in groups
        from halotools.utils import aggregation
        from scipy.stats import rankdata
        def function(x,key):
            return rankdata(-x[key],method='ordinal')-1
        mock = aggregation.add_members_property(mock, 'group_ID', 'rank', function, ['Mstar'])
    
        host = np.where(mock['rank']==0)[0]
        sub = np.where(mock['upid']!=0)[0]
        host_bool = (mock['rank']==0)
        sub_bool = (mock['rank']!=0)
    else:
        host = np.where(mock['upid']==-1)[0]
        sub = np.where(mock['upid']!=-1)[0]
        host_bool = (mock['upid']==-1)
        sub_bool = (mock['upid']!=-1)
    
    #star forming and quenched
    LHS = -11.0
    blue = np.where(mock['SSFR']>LHS)[0] #indices of blue galaxies
    red = np.where(mock['SSFR']<LHS)[0] #indicies of red galaxies

    #calculate the quenched fraction of centrals as a function of halo mass
    from f_prop import f_prop
    f_red_cen = f_prop(mass,bins,red,blue,host_bool)
    
    #calculate how many hosts are in each mass bin
    N = np.histogram(mass[host_bool],bins=bins)[0]
    
    return f_red_cen, N

def abundance_match(mock):

    #centrals and satellites
    host = np.where(mock['upid']==-1)[0]
    sub = np.where(mock['upid']!=-1)[0]
    host_bool = (mock['upid']==-1)
    sub_bool = (mock['upid']!=-1)
    N_sat = len(sub)
    N_cen = len(host)
    N_gal = len(host)+len(sub)
    
    #quick and dirty abundance matching
    Mam = np.zeros((len(host),))
    sorted_inds = np.argsort(mock['Mstar_all'][host])
    Mam[sorted_inds] = np.sort(mock['mvir'][host])
    Mam_all = np.zeros(len(mock))
    Mam_all[host] = Mam
    
    return Mam_all


if __name__ == '__main__':
    main()