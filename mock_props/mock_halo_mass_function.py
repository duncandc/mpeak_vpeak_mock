#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#fit the halo + sub halo mass function in Bolshoi.

#load packages
from __future__ import print_function, division
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import h5py
import sys

def main():

    #open halo catalogue
    Lbox=250.0
    filepath = cu.get_output_path() + 'processed_data/Multidark/Bolshoi/halo_catalogues/'
    halo_catalogue = 'hlist_1.00030.list'
    f =  h5py.File(filepath+halo_catalogue+'.hdf5', 'r')
    HC = f.get(halo_catalogue) #halo catalogue
    HC = np.array(HC)
    for name in HC.dtype.names: print(name)
    
    #make completeness cuts
    #require all (sub-)haloes to had had >100 particles at some point
    mp_bolshoi = 1.35e8
    mpeak_cut = mp_bolshoi*100.0
    keep = (HC['Mpeak']>mpeak_cut)
    HC = HC[keep]
    #keep = (HC['upid']==-1)
    #HC = HC[keep]
    
    #define halo property to use to make mass function
    halo_prop = 'Mpeak'
    
    #calculate cumulative mass function
    #N = np.cumsum(np.ones(len(HC)))
    #N = N/(Lbox**3.0)
    #M = np.sort(HC[halo_prop])[::-1]
    #M = np.log10(M)
    #plt.plot(M,N)
    #plt.yscale('log')
    #plt.show()
    
    #assign jackknife labels
    sample = np.vstack((HC['x'],HC['y'],HC['z'])).T
    Nsub=np.array([5.0,5.0,5.0])
    j_index, N_sub_vol = get_subvolume_labels(sample, Nsub, np.array([Lbox]*3))
    N_subs = get_subvolume_numbers(j_index,N_sub_vol)
    N_sub_vol = np.prod(Nsub).astype(np.int)
    
    #calculate mass function
    bins = np.arange(10.2,16.0,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    weights = np.array([1.0/(Lbox**3.0)]*len(HC))
    phi_full = mass_function(HC,halo_prop,bins,weights,use_log=True)
    phi_sub = np.zeros((N_sub_vol,len(phi_full)))
    for i in range(0,N_sub_vol):
        inds = (j_index==i+1)
        phi_sub[i,:] = mass_function(HC[inds], halo_prop, bins, weights[inds], use_log=True)
    
    cov = covariance_matrix(phi_sub,phi_full,N_sub_vol)
    err = np.sqrt(np.diagonal(cov))
    print(err)
    
    raw_counts = np.histogram(np.log10(HC[halo_prop]), bins=bins)[0]
    err_pos = np.sqrt(raw_counts)
    err_pos = (1.0/err_pos)*phi_full
    print(err_pos)
    
    print(err/err_pos)
    
    keep = (raw_counts>0)
    phi_full = phi_full[keep]
    err= err[keep]
    bin_centers = bin_centers[keep]
    
    #fit function to data
    from scipy.optimize import curve_fit
    log_err = 0.434*err/phi_full
    #single schechter funciton fit
    params, cov = curve_fit(Log_Schechter_function, bin_centers, np.log10(phi_full),\
                            p0=[0.012,12.4,-1.5, 1.0], sigma=log_err)
    dndm_halo_fit_1 = cu.schechter_function.Log_Super_Schechter(*params)
    print("params:", params)
    print("error:", np.sqrt(np.diagonal(cov)))
    red_chi_2 = cu.fitting.redchisqg(phi_full,dndm_halo_fit_1(bin_centers),deg=3,sd=err)
    print("Reduced chi squared value of fit: {0}".format(red_chi_2))
    
    """
    #double schechter funciton fit
    params, cov = curve_fit(Log_Double_Schechter_function, bin_centers, np.log10(phi),\
                            p0=[13.4,14.5,0.012,0.01,-1.5,-0.8], sigma=log_err)
    dndm_halo_fit_2 = cu.schechter_function.Log_Double_Schechter(*params)
    print(params)
    """
    """
    #special function fit
    params, cov = curve_fit(log_special_fit, bin_centers, np.log10(phi),\
                            p0=[0.0012,14.5,-2.0, 0.01, 12.95, 2.5], sigma=log_err)
    dndm_halo_fit_2 = lambda x: special_fit(x,*params)
    print(params)
    """
    
    
    #plot mass function
    m_sample=np.linspace(10,16,100)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    ax = fig.add_subplot(1,1,1)
    p1 = ax.errorbar(10**bin_centers,phi_full,yerr=err, fmt='.',color='black')
    p2, = ax.plot(10**m_sample,dndm_halo_fit_1(m_sample), color='black')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim([10**(-7),1])
    plt.xlabel(r'$M_{\rm peak}~[h^{-1}M_{\odot}]$')
    plt.ylabel(r'$\phi(M_{\rm peak})~[h^{3}{\rm Mpc}^{-3}\log(M_{\rm peak})^{-1}]$')
    plt.legend((p1,p2),("mock", "fit"), loc=3, fontsize=10, frameon=False, numpoints=1)
    plt.show()
    
    savepath = cu.get_plot_path()+'/analysis/central_quenching/'
    filename = 'halo_mass_function'
    fig.savefig(savepath+filename+'.pdf', dpi=250)
    
    
    """
    #plot difference in the mass function
    m_sample=np.linspace(10,16,100)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    ax = fig.add_subplot(1,1,1)
    ax.plot(bin_centers,dndm_halo_fit_2(bin_centers)/phi)
    ax.plot(bin_centers,dndm_halo_fit_1(bin_centers)/phi)
    plt.ylim([0.9,1.1])
    plt.show()
    """


def mass_function(HC,halo_prop,bins,weights,use_log=True):
    
    #calculate mass function
    bin_centers = (bins[:-1]+bins[1:])/2.0
    dm = (bins[1:]-bins[:-1])
    if use_log==True:
        counts = np.histogram(np.log10(HC[halo_prop]), bins=bins, weights=weights)[0]
    else:
        counts = np.histogram(HC[halo_prop], bins=bins, weights=weights)[0]
    phi = counts/dm
    
    return phi


def covariance_matrix(sub,full,N_sub_vol):
        """
        Calculate the covariance matrix.
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
        ind = np.random.randint(0,len(sample1),size=1)
       
        j_index_1 = inds[index_1[:,0],index_1[:,1],index_1[:,2]].astype(int)
        return j_index_1, N_sub_vol
    
def get_subvolume_numbers(j_index, N_sub_vol):
    temp = np.hstack((j_index,np.arange(1,N_sub_vol+1,1))) #kinda hacky, but need to every label to be in there at least once
    labels, N = np.unique(temp,return_counts=True)
    N = N-1 #remove the place holder I added two lines above.
    return N


def Log_Schechter_function(x, phi0, x0, alpha1, alpha2):
    x = np.asarray(x)
    x = x.astype(float)
    norm = np.log(10.0)*phi0
    val = norm*(10.0**((x-x0)*(1.0+alpha1)))*np.exp(-10.0**(alpha2*(x-x0)))
    return np.log10(val)

def log_special_fit(x, phi0, x0, alpha, a, b, c):
    x = np.asarray(x)
    x = x.astype(float)
    norm = np.log(10.0)*phi0
    val = norm*(10.0**((x-x0)*(1.0+alpha)))*np.exp(-10.0**(x-x0))
    
    fcor = a*np.exp(-1.0*(x-b)**2/(2.0*c**2))
    
    val = val*fcor
    
    return np.log10(val)

def special_fit(x, phi0, x0, alpha, a, b, c):
    x = np.asarray(x)
    x = x.astype(float)
    norm = np.log(10.0)*phi0
    val = norm*(10.0**((x-x0)*(1.0+alpha)))*np.exp(-10.0**(x-x0))
    
    fcor = a*np.exp(-1.0*(x-b)**2/(2.0*c**2))
    
    val = val*fcor
    
    return val

def Log_Double_Schechter_function(x, x1, x2, phi1, phi2, alpha1, alpha2):
    x = np.asarray(x)
    x = x.astype(float)
    norm = np.log(10.0)
    val = norm *\
        (np.exp(-10.0**(x - x1)) * 10.0**(x - x1) *\
            phi1 * (10.0**((x - x1) * alpha1)) +\
        np.exp(-10.0**(x - x2)) * 10.0**(x - x2) *\
            phi2 * (10.0**((x - x2) * alpha2)))
    return np.log10(val)


if __name__ == '__main__':
    main()