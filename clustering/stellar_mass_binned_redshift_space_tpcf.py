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
from halotools_devel.mock_observables import redshift_space_tpcf
from halotools_old.mock_observables import distant_observer
from numpy.lib.recfunctions import append_fields

def main():
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    
    if len(sys.argv)>1:
        catalogue = sys.argv[1]
        sm_low = float(sys.argv[2])
        sm_high = float(sys.argv[3])
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        sm_low=10.25
        sm_high=12.0
    
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
    mock = np.array(mock)
    print(mock.dtype.names)
    
    mock = append_fields(mock, 'redshift', np.zeros(len(mock)))
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    
    #from astropy.constants import c
    c = 299792.458 #speed of light in km/s
    from astropy import cosmology
    from scipy.interpolate import interp1d
    
    #get the peculiar velocity component along the line of sight direction
    v_los = mock['vz']
    
    #compute cosmological redshift
    y = np.linspace(0,10,1000)
    x = cosmo.comoving_distance(y).value
    f = interp1d(x, y, kind='cubic')
    finv = interp1d(y, x, kind='cubic')
    z_cos = f(mock['z'])
    
    #redshift is combination of cosmological and peculiar velocities
    z = z_cos+(v_los/c)*(1.0+z_cos)
    
    #reflect galaxies around redshift PBC
    flip = (z>f(period[2]))
    z[flip] = z[flip]-f(period[2])
    flip = (z<0)
    z[flip] = z[flip]+f(period[2])
    
    mock['redshift']=z
    
    mock['redshift'] = finv(mock['redshift'])
    print(mock['redshift'])
    
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
    rp_bins = np.linspace(-2.0,1.5,25)
    rp_bins = 10.0**rp_bins
    rp_bin_centers = (rp_bins[:-1]+rp_bins[1:])/2.0
    
    pi_bins = np.linspace(-2.0,1.5,25)
    pi_bins = 10.0**pi_bins
    pi_bin_centers = (pi_bins[:-1]+pi_bins[1:])/2.0
    
    #calculate correlation functions
    rp_bins = np.linspace(0,5,40)
    rp_bin_centers = (rp_bins[:-1]+rp_bins[1:])/2.0
    
    pi_bins = np.linspace(0,5,40)
    pi_bin_centers = (pi_bins[:-1]+pi_bins[1:])/2.0
    
    #full sample auto correlation
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection = (selection_1 & selection_2)
    sample1 = np.vstack((mock['x'],mock['y'],mock['redshift'])).T[selection]
    print("number of all galaxies: {0}".format(len(sample1)))
    
    result_all = redshift_space_tpcf(sample1, rp_bins, pi_bins,  period=period,
                                     do_auto=True, do_cross=False, estimator='Natural', 
                                     comm=None, max_sample_size=int(1e7))
    
    print(np.min(result_all),np.max(result_all))
    
    #rp_bin_centers = np.log10(rp_bin_centers)
    #pi_bin_centers = np.log10(pi_bin_centers)
    
    fig1, ax = plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True,figsize=(3.3, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.8, bottom=0.2, top=0.9)
    m = ax.imshow(np.flipud(np.log10(result_all).T),interpolation='nearest',aspect='auto',\
               extent=[rp_bin_centers[0],rp_bin_centers[-1],pi_bin_centers[0],pi_bin_centers[-1]])
    cax = fig1.add_axes([0.82, 0.2, 0.03, 0.7])
    fig1.colorbar(m, cax=cax)
    ax.set_xlabel('$r_p [Mpc/h]$')
    ax.set_ylabel('$\pi [Mpc/h]$')
    plt.show()
    
    fig1.savefig("/Users/duncan/Desktop/figure_1.pdf")

if __name__ == '__main__':
  main()