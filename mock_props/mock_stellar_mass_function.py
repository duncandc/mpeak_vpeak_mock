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
        catalogue = 'sm_9.5_s0.2_sfr_c-1.0_250'
    
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    
    print(mock.dtype.names)
    
    #measure stellar mass function
    bins = np.arange(9.5,12,0.1)
    dn, dm, err = raw_abundance(mock['Mstar'], 1.0/Lbox**3.0, bins, xlog=False)
    
    """
    #define fit to mass function, Baldry 2008.
    M0 = 10.648 + np.log10(0.7**2)
    phi1 = 4.26 * 10**(-3) / 0.7**3
    phi2 = 0.58 * 10**(-3) / 0.7**3
    alpha1 = -0.46
    alpha2 = -1.58
    dndm_gal = cu.schechter_function.Log_Double_Schechter(M0,M0,phi1,phi2,alpha1,alpha2)
    """
    #define stellar mass function
    #phi_1 = 7.97407978*10**(-3)
    #phi_2 = 1.10507295*10**(-3)
    #m_1 = 10.5040374
    #m_2 = 10.9083632
    #alpha_1 = -0.956108368
    #alpha_2 = -1.48157281
    phi_1 = 6.94559219e-03
    phi_2 = 2.76601643e-04
    m_1 = 1.05518723e+01
    m_2 = 1.08818830e+01
    alpha_1 = -1.20094175e+00
    alpha_2 = -7.40886755e-01
    dndm_gal = cu.schechter_function.Log_Double_Schechter(m_1, m_2, phi_1, phi_2, alpha_1, alpha_2)
    
    #plot stellar mass function
    msample=np.linspace(9,12,100)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    p1  = plt.errorbar(10**dm,dn,yerr=err,fmt='.',color='black')
    p2, = plt.plot(10**msample, dndm_gal(msample),color='black')
    plt.xlabel(r'$h^{-2}M_{*}~[M_{\odot}]$')
    plt.ylabel(r'$h^{3}\phi~[{\rm Mpc^{-3}}\log(M_{*})]$')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim([10**(-7),0.1])
    plt.legend((p1,p2),('mock','SDSS'), loc=3, fontsize=10, frameon=False, numpoints=1)
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'mock_stellar_mass_function'
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