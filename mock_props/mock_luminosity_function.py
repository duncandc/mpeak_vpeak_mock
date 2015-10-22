#Duncan Campbell
#March 2015
#Yale University
#Plot the stellar mass function from the mock

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    #open mock
    Lbox = 250.0
    if len(sys.argv)>1:
        catalogue = sys.argv[1]
    else:
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250_photo'
    
    #catalogue = 'Mr19_age_distribution_matching_mock'
    #filepath = '/Users/duncan/Documents/projects/output/processed_data/hearin_mocks/custom_catalogues/'
    
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    
    print(mock.dtype.names)
    
    #measure stellar mass function
    bins = np.arange(-24,-17,0.1)
    dn, dm, err = raw_abundance(mock['M_r,0.1'], 1.0/Lbox**3.0, bins, xlog=False)
    
    #luminosity function fit
    phi1=5.11595409e-03
    phi2=4.00805877e-03
    m1=-2.09667367e+01
    m2=-1.94722889e+01
    a1=-1.54378283e+00
    a2=6.61177120e-02
    dndm_gal = cu.schechter_function.Double_Mag_Schechter(phi1,phi2,m1,m2,a1,a2)
    
    #plot stellar mass function
    msample=np.arange(-24,-17,0.1)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    p1  = plt.errorbar(dm,dn,yerr=err,fmt='.',color='black')
    p2, = plt.plot(msample, dndm_gal(msample), '-',color='black')
    plt.xlabel(r'$M_{r,0.1}-5\log(h)$')
    plt.ylabel(r'$h^{3}\phi(M_{r,0.1}-5\log(h))~[{\rm Mpc^{-3}}{\rm dex}^{-1}]$')
    plt.yscale('log')
    #plt.xscale('log')
    plt.ylim([10**(-7),0.1])
    plt.xlim([-17,-24])
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'mock_luminosity_function'
    fig.savefig(savepath+filename+'.pdf', dpi=250)


def raw_abundance(x, weights, bins, xlog=True, monotonic_correction=False, show=False):
    """
    given objects with property 'x' and weights, return tabulated abundances.
    
    Parameters
    ==========
    x: array_like
    
    weights: array_like
    
    bins: array_like
    
    xlog: boolean, optional
        If True, use log values of x. 
    
    monotonic_correction: boolean, optional
        If True, attempt to correct for minor non-monotonicity
    
    show: boolean, optional
        plot abundance function 
    
    Returns
    =======
    dndx: dn, x, err
    """
    
    if xlog==True:
        x = np.log10(x)
    
    if np.shape(weights)==():
        weights = np.array([weights]*len(x))
    
    n = np.histogram(x, bins=bins, weights=weights)[0]
    bin_centers = (bins[:-1]+bins[1:])/2.0
    dx = bins[1:]-bins[:-1]
    dn = n/dx
    
    raw_counts = np.histogram(x, bins=bins)[0]
    
    #remove zeros
    keep = (raw_counts>0)
    dn = dn[keep]
    bin_centers = bin_centers[keep]
    raw_counts = raw_counts[keep]
    
    err = (1.0/np.sqrt(raw_counts))*dn
    
    if show==True:
        fig = plt.figure(figsize=(3.3,3.3))
        fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
        plt.errorbar(bin_centers, dn, yerr=err,fmt='.')
        plt.yscale('log')
        plt.show()
    
    return dn, bin_centers, err


if __name__ == '__main__':
    main()