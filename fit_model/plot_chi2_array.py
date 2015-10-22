#!/usr/bin/env python

#Duncan Campbell
#September 2015
#Yale University
#examine the covariance matrix of the correlation function for our fiducial mock

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.io import ascii
import custom_utilities as cu
from scipy.ndimage.filters import gaussian_filter
import sys

def main():

    N_top = 20

    filepath = './'
    filename = 'chinchilla_rho_3_sigma_2_jacknife_8_8_8_chi_2_array.npy'
    X = np.load(filepath+filename)
    
    X = X/(24.0*3.0*3.0)
    
    min_mask = (X==np.min(X))
    
    print(X[min_mask])
    
    inds = np.argsort(X, axis=None)
    inds = np.unravel_index(inds,(25,25))
    
    #mock parameters
    rhos = np.linspace(0.0,-1.0,25)
    sigmas = np.linspace(0.0,0.3,25)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0, wspace=0, left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.imshow(X.T, origin='lower',\
               extent=[sigmas[0],sigmas[-1],rhos[0],rhos[-1]],interpolation='nearest',\
               aspect='auto',cmap='hot', vmin=0.0, vmax=20.0)
    cbar = plt.colorbar()
    cbar.set_label(r'$\chi^2/{\rm dof}$')
    xx, yy = np.meshgrid(sigmas, rhos)
    CS = plt.contour(xx, yy, gaussian_filter(X.T,(1,1)), [3.0,4.0,6.0,8.0,10.0], colors='white')
    plt.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
    plt.plot([0.2],[-0.8],'x',color='red',ms=10)
    plt.plot([xx[min_mask.T]],[yy[min_mask.T]],'x',color='white',ms=10)
    plt.scatter(xx.T[inds[0][0:N_top],inds[1][0:N_top]], yy.T[inds[0][0:N_top],inds[1][0:N_top]],\
                lw=0, marker='.', color='grey')
    plt.xlabel(r'$\sigma_{\rm SMHM}$')
    plt.ylabel(r'$\rho_{\rm SSFR}$')
    plt.xlim([0,0.3])
    plt.ylim([0.0,-1.0])
    plt.show()
    
    fig.savefig('/Users/duncan/Desktop/chi_2_array.pdf', dpi=300)


if __name__ == '__main__':
    main()