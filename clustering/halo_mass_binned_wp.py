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
from halotools_old.mock_observables.two_point_functions import two_point_correlation_function

def main():
    
    N_threads = 1 #for calculating the TPCFs
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    
    if len(sys.argv)>1:
        catalogue = sys.argv[1]
        mvir_low = float(sys.argv[2])
        mvir_high = float(sys.argv[3])
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_cen_shuffle'
        #catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        mvir_low=12.0
        mvir_high=12.5
    
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
    rbins = np.linspace(-2.0,1.695,25)
    rbins = 10.0**rbins
    rbin_centers = (rbins[:-1]+rbins[1:])/2.0
    
    #full sample auto correlation
    selection_1 = (np.log10(mock['mvir'])>mvir_low)
    selection_2 = (np.log10(mock['mvir'])<mvir_high)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of all galaxies: {0}".format(len(sample1)))
    
    result_all = two_point_correlation_function(sample1, rbins, period=period,
                                   do_auto=True, do_cross=False, estimator='Natural', 
                                   N_threads=N_threads, comm=None,max_sample_size=int(1e7))
    
    #quenched sample auto correlation
    selection_1 = (np.log10(mock['mvir'])>mvir_low)
    selection_2 = (np.log10(mock['mvir'])<mvir_high)
    selection_1 = (selection_1 & selection_2)
    selection_2 = (red)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of red galaxies: {0}".format(len(sample1)))
    
    result_q = two_point_correlation_function(sample1, rbins, period=period,
                                   do_auto=True, do_cross=False, estimator='Natural', 
                                   N_threads=N_threads, comm=None,max_sample_size=int(1e7))
    
    #star-forming sample auto correlation
    selection_1 = (np.log10(mock['mvir'])>mvir_low)
    selection_2 = (np.log10(mock['mvir'])<mvir_high)
    selection_1 = (selection_1 & selection_2)
    selection_2 = (blue)
    selection = (selection_1 & selection_2 & host_bool)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    print("number of blue galaxies: {0}".format(len(sample1)))
    
    result_sf = two_point_correlation_function(sample1, rbins, period=period,
                                   do_auto=True, do_cross=False, estimator='Natural', 
                                   N_threads=N_threads, comm=None,max_sample_size=int(1e7))
    
    #save the correlation functions
    data_1 = Table([rbin_centers,result_all], names=['r', 'xi'])
    data_2 = Table([rbin_centers,result_q], names=['r', 'xi'])
    data_3 = Table([rbin_centers,result_sf], names=['r', 'xi'])
    
    filename_1 = catalogue+'_xi_all_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_1, savepath+filename_1)
    filename_2 = catalogue+'_xi_all_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_2, savepath+filename_2)
    filename_3 = catalogue+'_xi_all_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_3, savepath+filename_3)
    
    #calculate the projected correlation function
    from scipy import interpolate
    from scipy import integrate
    #use log(r) and log(xi) for interpolation
    x = np.log10(rbin_centers)
    y = np.log10(result_all+1)
    #xi_r_all = interpolate.interp1d(x, y)
    xi_r_all = interpolate.InterpolatedUnivariateSpline(x, y, k=1)
    y = np.log10(result_q+1)
    #xi_r_q = interpolate.interp1d(x, y)
    xi_r_q = interpolate.InterpolatedUnivariateSpline(x, y, k=1)
    y = np.log10(result_sf+1)
    #xi_r_sf = interpolate.interp1d(x, y)
    xi_r_sf = interpolate.InterpolatedUnivariateSpline(x, y, k=1)
    
    #define integrands
    def integrand_all(pi,rp):
        arg = np.sqrt(pi**2+rp**2)
        arg = np.log10(arg)
        return 10**xi_r_all(arg)-1
    def integrand_q(pi,rp):
        arg = np.sqrt(pi**2+rp**2)
        arg = np.log10(arg)
        return 10**xi_r_q(arg)-1
    def integrand_sf(pi,rp):
        arg = np.sqrt(pi**2+rp**2)
        arg = np.log10(arg)
        return 10**xi_r_sf(arg)-1
    
    #define projected radial bins
    rpbins = np.linspace(-1.0,np.log10(20.0),25)
    rpbins = 10.0**rpbins
    rpbin_centers = (rpbins[:-1]+rpbins[1:])/2.0
    #define arrays to hold results
    wp_q = np.zeros(len(rpbins)-1)
    wp_sf = np.zeros(len(rpbins)-1)
    wp_all = np.zeros(len(rpbins)-1)
    #define integration limits
    pi_min=0.0
    pi_max=40.0
    #do integrals
    for i in range(0,len(rpbins)-1):
        rp = rpbin_centers[i]
        print(i,rp,np.sqrt(rp**2+pi_max**2))
        wp_all[i] = 2.0*integrate.quad(integrand_all,pi_min,pi_max,args=(rp,))[0]
        wp_q[i] = 2.0*integrate.quad(integrand_q,pi_min,pi_max,args=(rp,))[0]
        wp_sf[i] = 2.0*integrate.quad(integrand_sf,pi_min,pi_max,args=(rp,))[0]
    
    #save projected correlation functions
    data_1 = Table([rpbin_centers,wp_all], names=['r', 'wp'])
    data_2 = Table([rpbin_centers,wp_q], names=['r', 'wp'])
    data_3 = Table([rpbin_centers,wp_sf], names=['r', 'wp'])
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    filename_1 = catalogue+'_wp_all_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_1, savepath+filename_1)
    filename_2 = catalogue+'_wp_q_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_2, savepath+filename_2)
    filename_3 = catalogue+'_wp_sf_halo_'+str(mvir_low)+'_'+str(mvir_high)+'.dat'
    ascii.write(data_3, savepath+filename_3)
    
    print(filename_1)

if __name__ == '__main__':
  main()