#Duncan Campbell
#March 2015
#Yale University
#shuffle mock to remove assembly bias, but preserve HOD

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    if len(sys.argv)>1:
        catalogue = sys.argv[1]
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        
    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(mock.dtype.names)

    hosts = (mock['upid']==-1)
    subs = (mock['upid']!=-1)
    
    bins = np.arange(10,15,0.1)
    shuffle_props = ['Mstar', 'SSFR']
    cen_shuffled_mock = shuffle_centrals(mock, hosts, bins, shuffle_props)

    bins = np.arange(10,15,0.1)
    shuffle_props = ['Mstar', 'SSFR']
    cen_sat_shuffled_mock = shuffle_satellites(cen_shuffled_mock, subs, bins, shuffle_props)
    
    bins = np.arange(10,15,0.1)
    shuffle_props = ['Mstar', 'SSFR']
    sat_shuffled_mock = shuffle_satellites(mock, subs, bins, shuffle_props)

    """
    selection = hosts
    plt.figure(figsize=(3.3,3.3))
    plt.scatter(cen_shuffled_mock['mvir'][selection],\
                cen_shuffled_mock['Mstar'][selection],\
                c=cen_shuffled_mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',\
                lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.xlim([10**10,10**15])
    plt.ylim([9.5,11.5])
    plt.xlabel(r'$M_{\rm vir}~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$M_{*}~[h^{-2}M_{\odot}]$')
    plt.show(block=True)
    
    selection = subs
    plt.figure(figsize=(3.3,3.3))
    plt.scatter(cen_sat_shuffled_mock['Mvir_host'][selection],\
                cen_sat_shuffled_mock['Mstar'][selection],\
                c=cen_sat_shuffled_mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',\
                lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.xlim([10**10,10**15])
    plt.ylim([9.5,11.5])
    plt.xlabel(r'$M_{\rm vir, host}~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$M_{*}~[h^{-2}M_{\odot}]$')
    plt.show(block=True)
    
    selection = hosts
    plt.figure(figsize=(3.3,3.3))
    plt.scatter(cen_sat_shuffled_mock['mvir'][selection],\
                cen_sat_shuffled_mock['Mstar'][selection],\
                c=cen_sat_shuffled_mock['SSFR'][selection],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',\
                lw=0, s=2, rasterized=True)
    plt.xscale('log')
    plt.xlim([10**10,10**15])
    plt.ylim([9.5,11.5])
    plt.xlabel(r'$M_{\rm vir}~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$M_{*}~[h^{-2}M_{\odot}]$')
    plt.show(block=True)
    """
    
    #check that central shuffle did not alter satellite props
    check = (cen_shuffled_mock[subs] == mock[subs])
    print("satellite props not changed in central shuffle:", np.all(check))
    
    #check that satellite shuffle did not alter central props
    check = (cen_sat_shuffled_mock[hosts] ==  cen_shuffled_mock[hosts])
    print("central props not changed in cenntral + satellite shuffle:", np.all(check))
    
     #check that satellite shuffle did not alter central props
    check = (sat_shuffled_mock[hosts] ==  mock[hosts])
    print("central props not changed in satellite shuffle:", np.all(check))
    
    #save mock
    savepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
    print('saving hdf5 version shuffled catalogue...')
    filename = catalogue+'_cen_shuffle'
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=cen_shuffled_mock)
    f.close()
    
    print('saving hdf5 version shuffled catalogue...')
    filename = catalogue+'_cen_sat_shuffle'
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=cen_sat_shuffled_mock)
    f.close()
    
    print('saving hdf5 version shuffled catalogue...')
    filename = catalogue+'_sat_shuffle'
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=sat_shuffled_mock)
    f.close()
    
    

def shuffle_centrals(mock, centrals, bins, shuffle_props, mock_prop='mvir', use_log=True):
    """
    Shuffle central galaxies amongst haloes.  
    Only one central per halo is allowed.
    """
    
    shuffled_mock = np.copy(mock)
    
    central_inds = np.where(centrals==True)[0]
    
    if use_log==True:
        inds = np.digitize(np.log10(mock[mock_prop][centrals]), bins=bins)
    else:
        inds = np.digitize(mock[mock_prop][centrals], bins=bins)
    
    
    for i in range(0,len(bins)-1):
        
        inds_in_bin = (inds==i+1)
        inds_in_bin = central_inds[inds_in_bin]
        
        shufled_inds_in_bin = np.random.permutation(inds_in_bin)
        
        for prop in shuffle_props:
            shuffled_mock[prop][shufled_inds_in_bin] = mock[prop][inds_in_bin]

    return shuffled_mock


def shuffle_satellites(mock, satellites, bins, shuffle_props, mock_prop_1='Mvir_host', mock_prop_2='Mpeak', use_log=True):
    """
    Shuffle satellite galaxies amongst haloes.
    
    First sort satellites by mock_prop_1. 
    Then shuffle satellite properties in bins of mock_prop_2.
    
    i.e. for satellites of the same host halo mass, shuffle satellites properties of the 
    same mpeak. 
    """
    
    shuffled_mock = np.copy(mock)
    
    #indicies of satellites
    satellite_inds = np.where(satellites==True)[0]
    
    #sort by host halo prop
    if use_log==True:
        inds = np.digitize(np.log10(mock[mock_prop_1][satellites]), bins=bins)
    else:
        inds = np.digitize(mock[mock_prop_1][satellites], bins=bins)
    
    #loop through bins of host halo prop
    for i in range(0,len(bins)-1):
        
        inds_in_bin = (inds==i+1)
        inds_in_bin = satellite_inds[inds_in_bin]
        
        #loop through bins of sub halo prop
        if len(inds_in_bin)>0:
            if use_log==True:
                sub_inds = np.digitize(np.log10(mock[mock_prop_2][inds_in_bin]), bins=bins)
            else:
                sub_inds = np.digitize(mock[mock_prop_2][inds_in_bin], bins=bins)
            for j in range(0,len(bins)-1):
                sub_inds_in_bin = (sub_inds==i+1)
                sub_inds_in_bin = inds_in_bin[sub_inds_in_bin]
            
                #shuffle satellite properties
                shufled_inds_in_bin = np.random.permutation(sub_inds_in_bin)
        
                for prop in shuffle_props:
                    shuffled_mock[prop][shufled_inds_in_bin] = mock[prop][sub_inds_in_bin]

    return shuffled_mock


if __name__ == '__main__':
    main()