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
from mpl_toolkits.mplot3d import Axes3D
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
    
    for i in range(284,360):
        print(i)
        fig = plt.figure(figsize=(10.0,10.0))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_axis_bgcolor('black')
        ax._axis3don = False
        #fig.subplots_adjust(left=0.2, right=0.82, bottom=0.2, top=0.9)
        r = [0, 250]
        from itertools import product, combinations
        for s, e in combinations(np.array(list(product(r,r,r))), 2):
            if np.sum(np.abs(s-e)) == r[1]-r[0]:
                ax.plot3D(*zip(s,e), color="white")
        p = ax.scatter(mock['x'],mock['y'],mock['z'],c=mock['SSFR'],\
                cmap='jet_r',vmax=-9.5, vmin=-12.5, marker='.',\
                lw=0, s=3, rasterized=False)
        ax.set_xlim([0,250])
        ax.set_ylim([0,250])
        ax.set_zlim([0,250])
        ax.view_init(azim=i)
    
        filepath = '/Users/duncan/Documents/projects/mpeak_vpeak_mock/movies/frames/movie_11/'
        filename = str(i).zfill(3)+'_mock_3D.png'
        plt.savefig(filepath+filename, dpi=400)
        plt.clf()
    

if __name__ == '__main__':
    main()