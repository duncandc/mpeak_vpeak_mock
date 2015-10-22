#!/usr/bin/env python

#Duncan Campbell
#May, 2015
#Yale University
#example script opening HDF5 mock

#The script can be run by typing "python example_open_mock.py" in the same directory as 
#    the mocks

#load packages
from __future__ import print_function
import numpy
import h5py


def main():

    filepath = './'
    catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
    
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open mock in 'read' mode
    mock = f.get(catalogue)
    
    #print out the column names in mock
    for i, name in enumerate(mock.dtype.names):
        print(i, name)
    
    #access the stellar mass column
    print(mock['Mstar'])
    
    #if you rather work with a numpy record array,
    mock = np.array(mock)


if __name__ == '__main__':
    main()