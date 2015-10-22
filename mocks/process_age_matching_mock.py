#!/usr/bin/env python

#Duncan Campbell
#April 8, 2014
#Yale University
#process the stellar mass age matching mock to be in the same format as my mocks

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    catalogue = 'sm9.8_age_matching_SFR_mock'
    sm_lim = 9.5
    
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/'
    print('opening mock catalogue:', catalogue+'.hdf5')
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    print(mock.dtype.names)

    new_mock = get_mock_array(mock)
    
    new_mock['id'] = mock['ID_halo']
    new_mock['upid'] = mock['ID_host']
    new_mock['x'] = mock['x']
    new_mock['y'] = mock['y']
    new_mock['z'] = mock['z']
    new_mock['vx'] = mock['Vx']
    new_mock['vy'] = mock['Vy']
    new_mock['vz'] = mock['Vz']
    new_mock['mvir'] = mock['M_vir']
    new_mock['vz'] = mock['Vz']
    new_mock['Mstar'] = mock['M_star']-(0.7**2)
    new_mock['SSFR'] = mock['SSFR']
    new_mock['Vpeak'] = mock['V_peak']
    
    #make completeness cut
    keep = (new_mock['Mstar']>sm_lim)
    new_mock = new_mock[keep]
    
    #save mock
    savepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
    print('saving hdf5 version of the extended catalogue...')
    filename = 'age_matching'
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=new_mock)
    f.close()


def get_mock_array(HC):
    #create structure to store mock
    dtype = [('id', '<i8'),('pid', '<i8'), ('upid', '<i8'),('mvir', '<f8'),\
             ('rvir', '<f8'),('vmax', '<f8'),('x', '<f8'), ('y', '<f8'), ('z', '<f8'),\
             ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),\
             ('Mvir_all', '<i8'),('M200b', '<i8'), ('M200c', '<i8'), ('M500c', '<i8'),('M2500c', '<i8'),\
             ('Macc', '<f8'),('Mpeak', '<f8'), ('Vacc', '<f8'), ('Vpeak', '<f8'),\
             ('First_Acc_Scale', '<f8'), ('Halfmass_Scale', '<f8'), ('Vmax@Mpeak', '<f8'),\
             ('Mstar','<f8'),('SSFR','<f8'),('Mstar_all','<f8'),('N_sat','<i8'),\
             ('Mvir_host','<f8'),('group_ID','<i8'),('SSFR_host','<f8'),('Rank','<i8'),\
             ('ejected', '<i8'),('abandoned', '<i8'),('ghost', '<i8'),\
             ('red', '<i8'), ('N', '<i8'), ('N_cen', '<f8')]
    dtype = np.dtype(dtype)
    
    #fill in some values
    mock = np.recarray((len(HC),), dtype=dtype)
    mock.fill(-99)
    
    return mock

if __name__ == '__main__':
    main()