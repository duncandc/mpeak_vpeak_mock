#!/usr/bin/env python

#Duncan Campbell
#January 29, 2015
#Yale University
#Calculate auto and cross correlation (SF and quenched) of stellar mass threshold samples 
#for central quenching mocks

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from halotools_old.mock_observables.two_point_functions import Delta_Sigma
from astropy.io import ascii
from astropy.table import Table

def main():
     
    if len(sys.argv)>1:
        catalogue = sys.argv[1]
        sm_bin_low = float(sys.argv[2])
        sm_bin_high = float(sys.argv[3])
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
        sm_bin_low = 10.0
        sm_bin_hgih = 10.5
    
    sm_bin = str(sm_bin_low)+'_'+str(sm_bin_high)
    
    #open galaxy mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    print('opening mock catalogue:', catalogue+'.hdf5')
    #open catalogue
    f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    
    #open dark matter particle data
    filepath_particles = cu.get_output_path() + 'processed_data/Multidark/Bolshoi/particle_catalogues/'
    part_catalogue = 'bolshoi_a1.0003_1e6_particles'
    f = h5py.File(filepath_particles+part_catalogue+'.hdf5', 'r') #open catalogue file
    particle_data = f.get(part_catalogue)
    
    particles = np.zeros((len(particle_data),3))
    particles[:,0]=particle_data['x']
    particles[:,1]=particle_data['y']
    particles[:,2]=particle_data['z']
    
    print("number of mater particles: {0}".format(len(particles)))
    
    #define star forming and quenched samples
    LHS = -11.0
    blue = (mock['SSFR']>LHS) #indices of blue galaxies
    red = (mock['SSFR']<LHS) #indicies of red galaxies
    
    #define radial bins
    rbins = np.linspace(-1,np.log10(2),15)
    rbins = 10.0**rbins
    rbins = np.insert(rbins,0,0.0)
    rbin_centers = (rbins[:-1]+rbins[1:])/2.0
    period = [250.0,250.0,250.0]
    
    #all galaxies
    selection = (mock['Mstar']>sm_bin_low) & (mock['Mstar']<sm_bin_high)
    centers = zip(mock['x'][selection],mock['y'][selection],mock['z'][selection])
    print("number of galaxies in selection: {0}".format(len(centers)))
    delta_sigma_all = Delta_Sigma(centers, particles, rbins, bounds=[-50.0,50.0], normal=[0.0,0.0,1.0],
                                  period=period, N_threads=1)
    
    #quenched galaxies
    selection_1 = (red)
    selection_2 = (mock['Mstar']>sm_bin_low) & (mock['Mstar']<sm_bin_high)
    selection = (selection_1 & selection_2)
    centers = zip(mock['x'][selection],mock['y'][selection],mock['z'][selection])
    print("number of galaxies in selection: {0}".format(len(centers)))
    delta_sigma_red = Delta_Sigma(centers, particles, rbins, bounds=[-50.0,50.0], normal=[0.0,0.0,1.0],
                                  period=period, N_threads=1)
    
    #star forming galaxies
    selection_1 = (blue)
    selection_2 = (mock['Mstar']>sm_bin_low) & (mock['Mstar']<sm_bin_high)
    selection = (selection_1 & selection_2)
    centers = zip(mock['x'][selection],mock['y'][selection],mock['z'][selection])
    print("number of galaxies in selection: {0}".format(len(centers)))
    delta_sigma_blue = Delta_Sigma(centers, particles, rbins, bounds=[-50.0,50.0], normal=[0.0,0.0,1.0],
                                  period=period, N_threads=1)
    
    #units?
    N_particles_tot = 2048.0**3.0
    N_particles_sub = len(particles)
    Mass_particles = 1.35*10.0**8.0
    Mass_particle_sub = N_particles_tot/N_particles_sub * Mass_particles
    #units #Msol/Mpc^2
    delta_sigma_blue = delta_sigma_blue*Mass_particle_sub
    delta_sigma_red = delta_sigma_red*Mass_particle_sub
    #now convert to Msol/pc^2
    delta_sigma_blue = delta_sigma_blue/((10.0**6.0)**2.0)
    delta_sigma_red = delta_sigma_red/((10.0**6.0)**2.0)
    
    
    data_1 = Table([rbin_centers,delta_sigma_all], names=['r', 'DeltaSigma'])
    data_2 = Table([rbin_centers,delta_sigma_red], names=['r', 'DeltaSigma'])
    data_3 = Table([rbin_centers,delta_sigma_blue], names=['r', 'DeltaSigma'])
    
    savepath = cu.get_output_path() + 'analysis/central_quenching/observables/'
    filename_1 = catalogue+'_DeltaSigma_sm_all_'+str(sm_bin)+'.dat'
    ascii.write(data_1, savepath+filename_1)
    filename_2 = catalogue+'_DeltaSigma_sm_q_'+str(sm_bin)+'.dat'
    ascii.write(data_2, savepath+filename_2)
    filename_3 = catalogue+'_DeltaSigma_sm_sf_'+str(sm_bin)+'.dat'
    ascii.write(data_3, savepath+filename_3)


if __name__ == '__main__':
    main()