#!/usr/bin/env python

#Duncan Campbell
#May, 2015
#Yale University
#add magnitudes and colors to stellar mass-ssfr mocks

#load packages
from __future__ import print_function
import numpy as np
import h5py
import custom_utilities as cu
import sys
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import append_fields

def main():

    if len(sys.argv)>1:
        catalogue = sys.argv[1]
        sm_lim = float(catalogue[3:6])
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        sm_lim = 9.5
    
    print(catalogue, sm_lim)
        
    #open mock
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    
    print(mock.dtype.names)
    
    mock_mr = np.zeros(len(mock))-99.0
    mock_color = np.zeros(len(mock))-99.0
    
    #open sdss catalogue
    filepath = cu.get_output_path() + 'processed_data/nyu_vagc/'
    filename = 'nyu_lss_mpa_vagc_dr7_sm_sample'
    f = h5py.File(filepath+filename+'.hdf5', 'r')
    dset = f.get(filename)
    dset = np.array(dset)
    
    #keep = (dset['Z']<0.04)
    #dset = dset[keep]
    
    #convert stellar mass measure to h=1 units
    dset['sm_MEDIAN'] = np.log10(10**dset['sm_MEDIAN']*0.7**2)
    
    #plt.figure()
    #plt.plot(dset['Z'],dset['sm_MEDIAN'],'.',ms=2)
    #plt.show()
    
    print(np.min(dset['sm_MEDIAN']),np.max(dset['sm_MEDIAN']))
    print(np.min(mock['Mstar']),np.max(mock['Mstar']))
    
    print(np.min(dset['sfr_MEDIAN']),np.max(dset['sfr_MEDIAN']))
    print(np.min(mock['SSFR']),np.max(mock['SSFR']))
    
    #define stellar mass bins
    sm_bins = np.arange(sm_lim,12.5,0.1)
    sdss_sm_binned = np.digitize(dset['sm_MEDIAN'],bins = sm_bins)
    mock_sm_binned = np.digitize(mock['Mstar'],bins = sm_bins)
    
    #define ssfr bins
    ssfr_bins = np.arange(-13.1,-8.0,0.2)
    sdss_ssfr_binned = np.digitize(dset['sfr_MEDIAN'],bins = ssfr_bins)
    mock_ssfr_binned = np.digitize(mock['SSFR'],bins = ssfr_bins)
    
    #loop through stellar mass bins
    for i,sm in enumerate(sm_bins):
        sdss_in_sm_bin = (sdss_sm_binned == i+1)
        mock_in_sm_bin = (mock_sm_binned == i+1)
        #loop through ssfr mass bins
        for j,ssfr in enumerate(ssfr_bins):
            sdss_in_ssfr_bin = (sdss_ssfr_binned == j+1)
            mock_in_ssfr_bin = (mock_ssfr_binned == j+1)
            
            mock_in_bin = (mock_in_ssfr_bin & mock_in_sm_bin)
            sdss_in_bin = (sdss_in_ssfr_bin & sdss_in_sm_bin)
            
            #number of objects in mock and sdss in bin
            N_mock = np.sum(mock_in_bin)
            N_sdss = np.sum(sdss_in_bin)
            
            if (N_mock>0) & (N_sdss>0):
                mags = dset['ABSMAG_r.nearest.model.z0.10'][sdss_in_bin]
                colors = (dset['ABSMAG_g.nearest.model.z0.10']-dset['ABSMAG_r.nearest.model.z0.10'])[sdss_in_bin]
                sample_inds = np.random.choice(np.arange(0,np.sum(sdss_in_bin),1),np.sum(mock_in_bin),replace=True) 
                
                mock_mr[mock_in_bin] = mags[sample_inds]
                mock_color[mock_in_bin] = colors[sample_inds]

            elif N_mock>0 & N_sdss==0:
                print(i,j,sm,ssfr,N_mock,N_sdss,N_mock,N_sdss)
                mock_mr[mock_in_bin] = -99
                mock_color[mock_in_bin] = -99
    
    #take care of galaxies which did not get properties assigned
    missing  = np.where(mock_mr==-99)[0]
    #find galaxies nearest to values
    for ind in missing:
        print(ind,mock['Mstar'][ind],mock['SSFR'][ind])
        d = np.fabs(dset['sm_MEDIAN']-mock['Mstar'][ind])+np.fabs(dset['sfr_MEDIAN']-mock['SSFR'][ind])
        sdss_ind = np.argmin(d)
        mock_mr[ind] = dset['ABSMAG_r.nearest.model.z0.10'][sdss_ind]
        mock_ssfr = (dset['ABSMAG_g.nearest.model.z0.10']-dset['ABSMAG_r.nearest.model.z0.10'])[sdss_ind]
        print(ind,mock['Mstar'][ind], mock['SSFR'][ind],\
              dset['sm_MEDIAN'][sdss_ind], dset['sfr_MEDIAN'][sdss_ind])
    
    mock = append_fields(mock,['M_r,0.1','(g-r)_0.1'],[mock_mr,mock_color])
    
    print('saving hdf5 version of the extended catalogue...')
    filename = catalogue+'_photo'
    savepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=mock)
    f.close()
    

if __name__ == '__main__':
    main()