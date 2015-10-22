#!/usr/bin/env python

#Author: Duncan Campbell
#January 28, 2015
#Yale University
#plot a slice of the mock box

#load packages
from __future__ import print_function
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    catalogue = 'sm_9.49_s0.2_sfr_c0.0_250'
    filepath_mock = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    
    #get slice in z
    show = (mock['z']<10.0)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.82, bottom=0.2, top=0.9)
    p = plt.scatter(mock['x'][show],mock['y'][show],c=mock['SSFR'][show],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',\
                lw=0, s=6, rasterized=False)
    plt.xlim([0,250])
    plt.ylim([0,250])
    plt.xlabel(r'$h^{-1} {\rm Mpc}$')
    plt.ylabel(r'$h^{-1} {\rm Mpc}$')
    cbar_ax = fig.add_axes([0.83, 0.2, 0.02, 0.7]) #xmin, ymin, +dx, +dy
    cbar = fig.colorbar(p, cax=cbar_ax)
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = catalogue+'_slice'
    fig.savefig(savepath+filename+'.png',dpi=400)
    

if __name__ == '__main__':
    main()