#!/usr/bin/env python

#Duncan Campbell
#January 29, 2015
#Yale University
#Calculate auto and cross correlation (SF and quenched) of stellar mass binned samples 
#for central quenching mocks

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from astropy.io import ascii
from astropy.table import Table
from halotools.mock_observables.clustering import wp
from halotools.mock_observables.pair_counters.rect_cuboid_pairs import xy_z_npairs

def main():
    
    N_threads = 4 #for calculating the TPCFs
    np.random.seed(0)
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    
    catalogue = 'sm_9.5_s0.2_sfr_c-0.8_Chinchilla_250'
    
    Nran = 10**6 #number of randoms...
    Nsub = np.array([8,8,8]) #jackknife sample splitting
    calculate_RR = False
    RR_filename = 'RR_8_8_8_10e6'
    
    if len(sys.argv)==3:
        sm_low=float(sys.argv[1])
        sm_high=float(sys.argv[2])
    else:
        sm_low=10.5
        sm_high=11.0
        
    print("Running with {} to {} stellar mass bin.".format(sm_low,sm_high))
    
    if catalogue[-3:]=='250':
        period = np.array([250.0,250.0,250.0])
    if catalogue[-3:]=='125':
        period = np.array([125.0,125.0,125.0])
    else: period = np.array([250.0,250.0,250.0])
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
    rp_bins = np.linspace(-1.0,1.4,25)
    rp_bins = 10.0**rp_bins
    rp_bin_centers = (rp_bins[:-1]+rp_bins[1:])/2.0
    
    pi_bins = np.linspace(0,50,50)
    pi_bin_centers = (pi_bins[:-1]+pi_bins[1:])/2.0
    
    
    #calculate jackknife samples
    labels = get_subvolume_labels(np.vstack((mock['x'],mock['y'],mock['z'])).T, Nsub,\
                                  np.array([250.0,250.0,250.0]))
    
    #calculate RR once!
    if calculate_RR==True:
        print('calculating RR')
        randoms = np.random.random((Nran,3))*250.0
        random_labels = get_subvolume_labels(randoms, Nsub, np.array([250.0,250.0,250.0]))
        r_selection = (random_labels!=1)
        RR = xy_z_npairs(randoms[r_selection], randoms[r_selection], rp_bins, pi_bins,\
                         period=period, N_threads=N_threads)
        RR = np.diff(np.diff(RR,axis=0),axis=1)
        NR = len(randoms[r_selection])
        np.save('./' + RR_filename,RR)
    else:
        print('loading RR')
        randoms = np.random.random((Nran,3))*250.0
        random_labels = get_subvolume_labels(randoms, Nsub, np.array([250.0,250.0,250.0]))
        r_selection = (random_labels!=1)
        NR = len(randoms[r_selection])
        RR = np.load('./' + RR_filename+'.npy')
    
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection = (selection_1 & selection_2)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    
    result_all = np.zeros((np.max(labels)+1,len(rp_bins)-1))
    result_all[0,:] = wp(sample1, rp_bins, pi_bins, period=period,
                         do_auto=True, do_cross=False, estimator='Natural',
                         N_threads=N_threads,max_sample_size=int(1e7))
    
    #loops over jackknife samples
    for i in range(1,np.max(labels)+1):
        print("all", i)
        #full sample auto correlation
        selection_3 = (labels!=i)
        selection = (selection_1 & selection_2 & selection_3)
        sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
        N1 = len(sample1)
        print("number of all galaxies: {0}".format(N1))
    
        DD = xy_z_npairs(sample1, sample1, rp_bins, pi_bins, period=period,\
                         N_threads=N_threads)
        DD = np.diff(np.diff(DD,axis=0),axis=1)
        xi = _TP_estimator(DD,None,RR,N1,N1,NR,NR,'Natural')
        result_all[i,:] = integrate_2D_xi(xi,pi_bins)
    
    
    result_q = np.zeros((np.max(labels)+1,len(rp_bins)-1))
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection_3 = (red)
    selection = (selection_1 & selection_2 & selection_3)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    result_q[0,:] = wp(sample1, rp_bins, pi_bins, period=period,
                       do_auto=True, do_cross=False, estimator='Natural', 
                       N_threads=N_threads,max_sample_size=int(1e7))
    
    for i in range(1,np.max(labels)+1):
        print("red", i)
        #quenched sample auto correlation
        selection_4 = (labels!=i)
        selection = (selection_1 & selection_2 & selection_3 & selection_4)
        sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
        N1 = len(sample1)
        print("number of red galaxies: {0}".format(N1))
    
        DD = xy_z_npairs(sample1, sample1, rp_bins, pi_bins, period=period,\
                         N_threads=N_threads)
        DD = np.diff(np.diff(DD,axis=0),axis=1)
        xi = _TP_estimator(DD,None,RR,N1,N1,NR,NR,'Natural')
        result_q[i,:] = integrate_2D_xi(xi,pi_bins)
    
    
    result_sf = np.zeros((np.max(labels)+1,len(rp_bins)-1))
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection_3 = (blue)
    selection = (selection_1 & selection_2 & selection_3)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
    result_sf[0,:] = wp(sample1, rp_bins, pi_bins, period=period,
                       do_auto=True, do_cross=False, estimator='Natural', 
                       N_threads=N_threads,max_sample_size=int(1e7))
    
    for i in range(1,np.max(labels)+1):
        print("blue", i)
        #star-forming sample auto correlation
        selection_4 = (labels!=i)
        selection = (selection_1 & selection_2 & selection_3 & selection_4)
        sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection]
        N1 = len(sample1)
        print("number of blue galaxies: {0}".format(N1))
    
        DD = xy_z_npairs(sample1, sample1, rp_bins, pi_bins, period=period,\
                         N_threads=N_threads)
        DD = np.diff(np.diff(DD,axis=0),axis=1)
        xi = _TP_estimator(DD,None,RR,N1,N1,NR,NR,'Natural')
        result_sf[i,:] = integrate_2D_xi(xi,pi_bins)
    
    result_cross = np.zeros((np.max(labels)+1,len(rp_bins)-1))
    selection_1 = (mock['Mstar']>sm_low)
    selection_2 = (mock['Mstar']<sm_high)
    selection_3 = (blue)
    selection_a = (selection_1 & selection_2 & selection_3)
    sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection_a]
    selection_3 = (red)
    selection_b = (selection_1 & selection_2 & selection_3)
    sample2 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection_b]
    result_cross[0,:] = wp(sample1, rp_bins, pi_bins, sample2=sample2, period=period,
                           do_auto=False, do_cross=True, estimator='Natural',
                           N_threads=N_threads, max_sample_size=int(1e7))
    
    for i in range(1,np.max(labels)+1):
        print("cross", i)
        #star-forming sample auto correlation
        selection_4 = (labels!=i)
        selection_aa = (selection_a & selection_4)
        selection_bb = (selection_b & selection_4)
        sample1 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection_aa]
        sample2 = np.vstack((mock['x'],mock['y'],mock['z'])).T[selection_bb]
        N1 = len(sample1)
        N2 = len(sample2)
        print("number of blue by red galaxies: {0} x {1}".format(N1, N2))
    
        DD = xy_z_npairs(sample1, sample2, rp_bins, pi_bins, period=period,\
                         N_threads=N_threads)
        DD = np.diff(np.diff(DD,axis=0),axis=1)
        xi = _TP_estimator(DD,None,RR,N1,N2,NR,NR,'Natural')
        result_cross[i,:] = integrate_2D_xi(xi,pi_bins)
    
    all_full = result_all[0,:]
    all_sub = result_all[1:,:]
    
    q_full = result_q[0,:]
    q_sub = result_q[1:,:]
    
    sf_full = result_sf[0,:]
    sf_sub = result_sf[1:,:]
    
    cross_full = result_cross[0,:]
    cross_sub = result_cross[1:,:]
    
    #save raw results
    filename_1 = catalogue+'_wp_fiducial_all_'+str(sm_low)+'_'+str(sm_high)
    np.save(filename_1+'_result', result_all)
    filename_2 = catalogue+'_wp_fiducial_q_'+str(sm_low)+'_'+str(sm_high)
    np.save(filename_2+'_result', result_q)
    filename_3 = catalogue+'_wp_fiducial_sf_'+str(sm_low)+'_'+str(sm_high)
    np.save(filename_3+'_result', result_sf)
    filename_4 = catalogue+'_wp_fiducial_cross_'+str(sm_low)+'_'+str(sm_high)
    np.save(filename_4+'_result', result_cross)
    
    #calculate covariance matrices
    cov_all = covariance_matrix(all_sub,all_full,np.max(labels))
    cov_red = covariance_matrix(q_sub,q_full,np.max(labels))
    cov_blue = covariance_matrix(sf_sub,sf_full,np.max(labels))
    cov_cross = covariance_matrix(cross_sub,cross_full,np.max(labels))
    
    #save projected correlation functions
    data_1 = Table([rp_bin_centers,all_full], names=['r', 'wp'])
    data_2 = Table([rp_bin_centers,q_full], names=['r', 'wp'])
    data_3 = Table([rp_bin_centers,sf_full], names=['r', 'wp'])
    data_4 = Table([rp_bin_centers,cross_full], names=['r', 'wp'])
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    
    filename_1 = catalogue+'_wp_fiducial_all_'+str(sm_low)+'_'+str(sm_high)
    ascii.write(data_1, savepath+filename_1+'.dat')
    np.save(filename_1+'_cov', cov_all)
    
    filename_2 = catalogue+'_wp_fiducial_q_'+str(sm_low)+'_'+str(sm_high)
    ascii.write(data_2, savepath+filename_2+'.dat')
    np.save(filename_2+'_cov', cov_red)
    
    filename_3 = catalogue+'_wp_fiducial_sf_'+str(sm_low)+'_'+str(sm_high)
    ascii.write(data_3, savepath+filename_3+'.dat')
    np.save(filename_3+'_cov', cov_blue)
    
    filename_4 = catalogue+'_wp_fiducial_cross_'+str(sm_low)+'_'+str(sm_high)
    ascii.write(data_4, savepath+filename_4+'.dat')
    np.save(filename_4+'_cov', cov_cross)
    
    plt.figure()
    plt.errorbar(rp_bin_centers, all_full, yerr=np.sqrt(np.diagonal(cov_all)),\
                 color='black')
    plt.errorbar(rp_bin_centers, q_full, yerr=np.sqrt(np.diagonal(cov_red)),\
                 color='red')
    plt.errorbar(rp_bin_centers, sf_full, yerr=np.sqrt(np.diagonal(cov_blue)),\
                 color='blue')
    plt.errorbar(rp_bin_centers, cross_full, yerr=np.sqrt(np.diagonal(cov_cross)),\
                 color='green')
    plt.xscale('log')
    plt.yscale('log')
    plt.show(block=False)


def get_subvolume_labels(sample1, Nsub, Lbox):
        """
        Split volume into subvolumes, and tag points in subvolumes with integer labels for 
        use in the jackknife calculation.
        
        note: '0' tag should be reserved. It is used in the jackknife calculation to mean
        'full sample'
        """
        
        dL = Lbox/Nsub # length of subvolumes along each dimension
        N_sub_vol = np.prod(Nsub) # total the number of subvolumes
        inds = np.arange(1,N_sub_vol+1).reshape(Nsub[0],Nsub[1],Nsub[2])
    
        #tag each particle with an integer indicating which jackknife subvolume it is in
        #subvolume indices for the sample1 particle's positions
        index_1 = np.floor(sample1/dL).astype(int)
        j_index_1 = inds[index_1[:,0],index_1[:,1],index_1[:,2]].astype(int)
        
        return j_index_1

def covariance_matrix(sub,full,N_sub_vol):
        """
        Calculate the full covariance matrix.
        """
        Nr = full.shape[0] # Nr is the number of radial bins
        cov = np.zeros((Nr,Nr)) # 2D array that keeps the covariance matrix 
        after_subtraction = sub - np.mean(sub,axis=0)
        tmp = 0
        for i in range(Nr):
            for j in range(Nr):
                tmp = 0.0
                for k in range(N_sub_vol):
                    tmp = tmp + after_subtraction[k,i]*after_subtraction[k,j]
                cov[i,j] = (((N_sub_vol-1)/N_sub_vol)*tmp)
    
        return cov

#integrate the redshift space TPCF to get w_p
def integrate_2D_xi(x,pi_bins):
    return 2.0*np.sum(x*np.diff(pi_bins),axis=1)

def _TP_estimator(DD,DR,RR,ND1,ND2,NR1,NR2,estimator):
    """
    private internal function.
    
    two point correlation function estimator
    
    note: jackknife_tpcf uses its own intenral version, this is not totally ideal.
    """
    if estimator == 'Natural':
        factor = ND1*ND2/(NR1*NR2)
        #DD/RR-1
        xi = (1.0/factor)*DD/RR - 1.0
    elif estimator == 'Davis-Peebles':
        factor = ND1*ND2/(ND1*NR2)
        #DD/DR-1
        xi = (1.0/factor)*DD/DR - 1.0
    elif estimator == 'Hewett':
        factor1 = ND1*ND2/(NR1*NR2)
        factor2 = ND1*NR2/(NR1*NR2)
        #(DD-DR)/RR
        xi = (1.0/factor1)*DD/RR - (1.0/factor2)*DR/RR
    elif estimator == 'Hamilton':
        #DDRR/DRDR-1
        xi = (DD*RR)/(DR*DR) - 1.0
    elif estimator == 'Landy-Szalay':
        factor1 = ND1*ND2/(NR1*NR2)
        factor2 = ND1*NR2/(NR1*NR2)
        #(DD - 2.0*DR + RR)/RR
        xi = (1.0/factor1)*DD/RR - (1.0/factor2)*2.0*DR/RR + 1.0
    else: 
        raise ValueError("unsupported estimator!")
    return xi


if __name__ == '__main__':
  main()