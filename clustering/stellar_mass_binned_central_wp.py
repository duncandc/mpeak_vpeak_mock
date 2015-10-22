#!/usr/bin/env python

#Duncan Campbell
#January 29, 2015
#Yale University
#Calculate auto and cross correlation (SF and quenched) of stellar mass binned samples 
#for central quenching mocks

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from astropy.io import ascii
from astropy.table import Table
from halotools.mock_observables.clustering import wp

def main():
    
    N_threads = 1 #for calculating the TPCFs
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    
    if len(sys.argv)>1:
        catalogue = sys.argv[1]
        sm_low = float(sys.argv[2])
        sm_high = float(sys.argv[3])
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        sm_low=9.5
        sm_high=9.8
    
    if catalogue[-3:]=='250':
        period = [250.0,250.0,250.0]
    if catalogue[-3:]=='125':
        period = [125.0,125.0,125.0]
    else: period = [250.0,250.0,250.0]
    print(period)
    
    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    print('opening mock catalogue:', catalogue+'.hdf5')
    #open catalogue
    f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    print(mock.dtype.names)
    
    #centrals and satellites
    host = np.where(mock['upid']==-1)[0]
    sub = np.where(mock['upid']!=-1)[0]
    host_bool = (mock['upid']==-1)
    sub_bool = (mock['upid']!=-1)
    N_sat = len(sub)
    N_cen = len(host)
    N_gal = len(host)+len(sub)
    print(N_gal, N_cen, N_sat)
    
    #star forming and quenched
    LHS = -11.0
    blue  = (mock['SSFR']>LHS) #indices of blue galaxies
    red   = (mock['SSFR']<LHS) #indicies of red galaxies
    
    #calculate correlation functions
    rp_bins = np.linspace(-2.0,1.4,25)
    rp_bins = 10.0**rp_bins
    rp_bin_centers = (rp_bins[:-1]+rp_bins[1:])/2.0
    
    pi_bins = np.linspace(0,50,50)
    pi_bin_centers = (pi_bins[:-1]+pi_bins[1:])/2.0
    
    #full sample auto correlation
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of all galaxies: {0}".format(len(sample1)))
    
    result_all = wp(sample1, rp_bins, pi_bins, period=period,
                    do_auto=True, do_cross=False, estimator='Natural', 
                    N_threads=N_threads, max_sample_size=int(1e7))
    
    #quenched sample auto correlation
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection_1 = (selection_1 & selection_2)
    selection_2 = (red)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of red galaxies: {0}".format(len(sample1)))
    
    result_q = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads, max_sample_size=int(1e7))
    
    #star-forming sample auto correlation
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection_1 = (selection_1 & selection_2)
    selection_2 = (blue)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of blue galaxies: {0}".format(len(sample1)))
    
    result_sf = wp(sample1, rp_bins, pi_bins, period=period,
                  do_auto=True, do_cross=False, estimator='Natural', 
                  N_threads=N_threads, max_sample_size=int(1e7))
    
    #save projected correlation functions
    data_1 = Table([rp_bin_centers,result_all], names=['r', 'wp'])
    data_2 = Table([rp_bin_centers,result_q], names=['r', 'wp'])
    data_3 = Table([rp_bin_centers,result_sf], names=['r', 'wp'])
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    filename_1 = catalogue+'_wp_central_all_'+str(sm_low)+'_'+str(sm_high)+'.dat'
    ascii.write(data_1, savepath+filename_1)
    filename_2 = catalogue+'_wp_central_q_'+str(sm_low)+'_'+str(sm_high)+'.dat'
    ascii.write(data_2, savepath+filename_2)
    filename_3 = catalogue+'_wp_central_sf_'+str(sm_low)+'_'+str(sm_high)+'.dat'
    ascii.write(data_3, savepath+filename_3)
    

if __name__ == '__main__':
  main()