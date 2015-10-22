#!/usr/bin/env python

#Duncan Campbell
#March, 2015
#Yale University
#make a SHAM mock with stellar masses and SSFRs.

#load packages
from __future__ import print_function, division
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import sys
import h5py
from astropy.cosmology import FlatLambdaCDM
import SHAM
from astropy.io import ascii
from astropy import table
from scipy import interpolate
import time
from abundance_matching.AM import AM_nonparam_1
from abundance_matching import SMHM_model
from abundance_matching.make_mocks import make_SHAM_mock
from rank_correlation_matching.rank_correlate import correlate

def main():

    #set random seed to make mocks reproducible
    np.random.seed(1)

    #get parameters to make mock
    if len(sys.argv)<5:
        dm_sim = 'Bolshoi_250'
        sfr_cor = -1.0
        sm_lim = 9.5
        sm_scat = 0.2
    else:
        sfr_cor = float(sys.argv[1])
        sm_lim = float(sys.argv[2])
        sm_scat = float(sys.argv[3])
        dm_sim = sys.argv[4]
    
    if dm_sim =='Bolshoi_250': box_size=250.0
    if dm_sim == 'Chinchilla_250': box_size=250.0
    elif dm_sim == 'Chinchilla_125': box_size=125.0
    
    mock = return_mock(sm_lim = sm_lim, sigma = sm_scat, rho = sfr_cor, sim = dm_sim)
    
    #save mock
    savepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    
    print('saving hdf5 version of the extended catalogue...')
    filename = 'sm_'+str(sm_lim)+'_s'+str(sm_scat)+'_sfr_c'+str(sfr_cor)+'_'+dm_sim
    print(filename)
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=mock)
    f.close()


def return_mock(sm_lim=9.5, sigma=0.2, rho=-0.1, sim='Bolshoi'):
    
    #get parameters to make mock
    dm_sim = sim
    sfr_cor = rho
    sm_lim = sm_lim
    sm_scat = sigma
    
    sm_mock_key = 'Mpeak'
    ssfr_mock_key = 'Vpeak'
    
    #load halo catalogue
    HC, Lbox = get_mock(dm_sim)
    box_size=Lbox
    vol = Lbox**3.0
    
    #load sdss catalogue
    GC = get_sdss_sample()
    
    #get mock data structure
    mock = get_mock_array(HC)
    
    #define central moment SMHM distribution
    from scipy.stats import norm
    def sigma2_xy(y):
        return sm_scat
    def P_x(y,mu_xy,sigma2_xy=sigma2_xy):
        mu_x = mu_xy(y)
        sigma2_x = sigma2_xy(y)
        p = norm(loc=mu_x, scale=sigma2_x)
        return p
    
    #define halo mass function
    phi_star = 9.95231920*10**(-5)
    m_star = 13.9833298
    alpha1 = -1.85781918
    alpha2 = 0.438785492
    dndm_halo = cu.schechter_function.Log_Super_Schechter(phi_star, m_star, alpha1, alpha2)
    
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
    
    #solve for the mean SMHM relation
    dm_star = np.linspace(5,12.5,100)
    dm_halo = np.linspace(10,15.5,100)
    dn_star = dndm_gal(dm_star)
    dn_halo = dndm_halo(dm_halo)
    
    if sm_scat==0.0:
        SMHM = AM_nonparam_1(dn_star, dm_star, dn_halo, dm_halo,\
                         None, y_min = np.amin(dm_halo),\
                         y_max = np.amax(dm_halo), ny=30, tol=0.001)
    else:    
        SMHM = AM_nonparam_1(dn_star, dm_star, dn_halo, dm_halo,\
                         P_x, y_min = np.amin(dm_halo),\
                         y_max = np.amax(dm_halo), ny=30, tol=0.001)
    
    #apply the mean SMHM relation
    P = lambda y: P_x(y, mu_xy=SMHM)
    
    #trim mock
    keep  = (np.log10(mock['Mpeak'])>10.0)
    mock = mock[keep]
    
    #make mock
    mock = make_SHAM_mock(mock, P, mock_prop=sm_mock_key, gal_prop='Mstar',\
                          use_log_mock_prop=True)
    
    #make completeness cut
    keep = (mock['Mstar']>sm_lim)
    mock = mock[keep]
    
    #add SSFRs
    bins = np.arange(sm_lim,12.5,0.1)
    mock = add_ssfrs(mock, GC, bins)
    
    #correlate SSFR and Vpeak
    mock = correlate(mock, 'Mstar', ssfr_mock_key, 'SSFR', bins=bins, rho=sfr_cor)
    
    #identify host and sub haloes
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    
    #star forming and quenched
    LHS = -11.0
    blue = (mock['SSFR']>LHS) #indices of blue galaxies
    red = (mock['SSFR']<=LHS) #indicies of red galaxies
    mock['red'] = 0
    mock['red'][red]=1
    
    #identify ejected haloes
    if np.all(mock['First_Acc_Scale']==-99):
        mock['ejected']=-99
    else:
        ejected_haloes = (mock['upid']==-1) & (mock['First_Acc_Scale']<1.0)
        mock['ejected'] = 0
        mock['ejected'][ejected_haloes] = 1
    
    #identify abandoned satellites and satellites with ghost central
    host_ids = np.unique(mock['id'][host])
    has_host = np.in1d(mock['upid'], host_ids)
    had_host = np.in1d(mock['upid'], mock['id'])
    has_no_host = (has_host==False)
    had_no_host = (had_host==False)
    
    abandoned = (has_no_host & had_host)
    ghost = (has_no_host & had_no_host & sub)
    
    mock['abandoned'] = 0
    mock['abandoned'][abandoned] = 1
    mock['ghost'] = 0
    mock['ghost'][ghost] = 1
    
    ####add some group properties#########################################################
    #groupID
    mock['group_ID'] = mock['upid']
    mock['group_ID'][host] = mock['id'][host]
    
    """
    from halotools.utils import aggregation
    def broadcast(x):
        central = (x['upid']==-1)
        return x['mvir'][central]
    mock = aggregation.add_members_property(mock, 'group_ID', 'Mvir_host', broadcast)
    
    def broadcast(x):
        central = (x['upid']==-1)
        return x['SSFR'][central]
    mock = aggregation.add_members_property(mock, 'group_ID', 'SSFR_host', broadcast)
    
    def logsum(x):
        vals = 10.0**(x['Mstar'])
        return np.log10(np.sum(vals))
    mock = aggregation.add_members_property(mock, 'group_ID', 'Mstar_all', logsum)
    
    def N(x):
        return len(x)
    mock = aggregation.add_members_property(mock, 'group_ID', 'N', N)
    
    def N(x):
        cen = (x['upid']==-1)
        return np.sum(cen)
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_cen', N)
    
    def N(x):
        sats = (x['upid']!=-1)
        return np.sum(sats)
    mock = aggregation.add_members_property(mock, 'group_ID', 'N_sat', N)
    
    from scipy.stats import rankdata
    def rank(x):
        return rankdata(-1.0*x['Mstar'])-1
    mock = aggregation.add_members_property(mock, 'group_ID', 'Rank', rank)
    ######################################################################################
    """
    
    return mock


def return_base_mock(sm_lim=9.5, sigma=0.2, sim='Bolshoi'):
    
    #get parameters to make mock
    dm_sim = sim
    sm_lim = sm_lim
    sm_scat = sigma
    
    sm_mock_key = 'Mpeak'
    ssfr_mock_key = 'Vpeak'
    
    #load halo catalogue
    HC, Lbox = get_mock(dm_sim)
    box_size=Lbox
    vol = Lbox**3.0
    
    #load sdss catalogue
    GC = get_sdss_sample()
    
    #get mock data structure
    mock = get_mock_array(HC)
    
    #define central moment SMHM distribution
    from scipy.stats import norm
    def sigma2_xy(y):
        return sm_scat
    def P_x(y,mu_xy,sigma2_xy=sigma2_xy):
        mu_x = mu_xy(y)
        sigma2_x = sigma2_xy(y)
        p = norm(loc=mu_x, scale=sigma2_x)
        return p
    
    #define halo mass function
    phi_star = 9.95231920*10**(-5)
    m_star = 13.9833298
    alpha1 = -1.85781918
    alpha2 = 0.438785492
    dndm_halo = cu.schechter_function.Log_Super_Schechter(phi_star, m_star, alpha1, alpha2)
    
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
    
    #solve for the mean SMHM relation
    dm_star = np.linspace(5,12.5,100)
    dm_halo = np.linspace(10,15.5,100)
    dn_star = dndm_gal(dm_star)
    dn_halo = dndm_halo(dm_halo)
    
    if sm_scat==0.0:
        SMHM = AM_nonparam_1(dn_star, dm_star, dn_halo, dm_halo,\
                         None, y_min = np.amin(dm_halo),\
                         y_max = np.amax(dm_halo), ny=30, tol=0.001)
    else:    
        SMHM = AM_nonparam_1(dn_star, dm_star, dn_halo, dm_halo,\
                         P_x, y_min = np.amin(dm_halo),\
                         y_max = np.amax(dm_halo), ny=30, tol=0.001)
    
    #apply the mean SMHM relation
    P = lambda y: P_x(y, mu_xy=SMHM)
    
    #trim mock
    keep  = (np.log10(mock['Mpeak'])>10.0)
    mock = mock[keep]
    
    #make mock
    mock = make_SHAM_mock(mock, P, mock_prop=sm_mock_key, gal_prop='Mstar',\
                          use_log_mock_prop=True)
    
    #make completeness cut
    keep = (mock['Mstar']>sm_lim)
    mock = mock[keep]
    
    #add SSFRs
    bins = np.arange(sm_lim,12.5,0.1)
    mock = add_ssfrs(mock, GC, bins)
    
    return mock


def return_quick_ssfr_mock(mock, rho=-0.1):
    
    #get parameters to make mock
    sfr_cor = rho
    
    #correlate SSFR and Vpeak
    bins = np.arange(np.min(mock['Mstar']),12.5,0.1)
    mock = correlate(mock, 'Mstar', 'Vpeak', 'SSFR', bins=bins, rho=sfr_cor)
    
    ####add some group properties#########################################################
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    #groupID
    mock['group_ID'] = mock['upid']
    mock['group_ID'][host] = mock['id'][host]
    
    from halotools.utils import aggregation
    def logsum(x):
        vals = 10.0**(x['Mstar'])
        return np.log10(np.sum(vals))
    mock = aggregation.add_members_property(mock, 'group_ID', 'Mstar_all', logsum)
    
    return mock


def emperical_abundance(mock, halo_prop, bins, Lbox = 250.0, use_log = True):
    
    bins = np.arange(10.2,16.0,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    dm = (bins[1:]-bins[:-1])
    weights = [1.0/(Lbox**3.0)]*len(mock)
    if use_log==True:
        counts = np.histogram(np.log10(mock[halo_prop]), bins=bins, weights=weights)[0]
    else:
        counts = np.histogram(mock[halo_prop], bins=bins, weights=weights)[0]
    phi = counts/dm
    
    #remove 0 values
    keep = counts>0
    phi = phi[keep]
    bin_centers = bin_centers[keep]
    
    return phi, bin_centers


def add_ssfrs(mock,GC,bins):
    """
    Add SSFRs to the mock
    """
    
    sm_key = 'sm_MEDIAN'
    Mstar_sdss = GC[sm_key]
    Mstar_mock = mock['Mstar']
    
    #now we need galaxy colors
    sdss_sfr = GC['sfr_MEDIAN']
    mock_sfr = np.ones(len(mock))*-99.0
    
    sm_bins = bins
    sm_binned_sdss_gals = np.digitize(Mstar_sdss,bins=sm_bins)
    sm_binned_mock_gals = np.digitize(Mstar_mock,bins=sm_bins)
    
    from scipy.interpolate import interp1d
    for i in range(0,len(sm_bins)-1):
        #print(i,sm_bins[i],sm_bins[i+1])
        
        mock_inds = np.where(sm_binned_mock_gals==i+1)[0]
        sdss_inds = np.where(sm_binned_sdss_gals==i+1)[0]
        
        if len(mock_inds)==0:
            continue
        
        N_mock = len(mock_inds)
        N_sdss = len(sdss_inds)
       
        sub_sfr_colors = sdss_sfr[sdss_inds]
        sorted_sfr = np.sort(sub_sfr_colors)
        
        count = np.ones(N_sdss)
        counts = np.cumsum(count)-1
        
        if N_sdss<1:
            mock_sfr[mock_inds]=-99.9
        else:
            normalized_counts = counts/float(counts[-1])
            f = interp1d(normalized_counts,sorted_sfr,kind='linear')
       
            sample = np.random.random(N_mock)
       
            sample_sfr = f(sample)
        
            mock_sfr[mock_inds] = sample_sfr
    
    mock['SSFR'] =  mock_sfr
    
    return mock


def get_mock_array(HC):
   #create structure to store mock
    dtype = [('id', '<i8'),('pid', '<i8'), ('upid', '<i8'),('mvir', '<f8'),\
             ('rvir', '<f8'),('vmax', '<f8'),('x', '<f8'), ('y', '<f8'), ('z', '<f8'),\
             ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'),\
             ('Mvir_all', '<i8'),('M200b', '<i8'), ('M200c', '<i8'), ('M500c', '<i8'),('M2500c', '<i8'),\
             ('Macc', '<f8'),('Mpeak', '<f8'), ('Vacc', '<f8'), ('Vpeak', '<f8'),\
             ('First_Acc_Scale', '<f8'), ('Halfmass_Scale', '<f8'), ('Vmax@Mpeak', '<f8'),\
             ('Mstar','<f8'),('SSFR','<f8'),('Mstar_all','<f8'),\
             ('Mvir_host','<f8'),('group_ID','<i8'),('SSFR_host','<f8'),('Rank','<i8'),\
             ('ejected', '<i8'),('abandoned', '<i8'),('ghost', '<i8'),\
             ('red', '<i8'), ('N', '<i8'), ('N_cen', '<f8'), ('N_sat', '<f8')]
    dtype = np.dtype(dtype)
    
    #fill in some values
    mock = np.recarray((len(HC),), dtype=dtype)
    mock.fill(-99)
    for name in mock.dtype.names[:-15]:
        if name not in HC.dtype.names:
            mock[name]=-99
        else:
            mock[name] = HC[name]
    
    return mock


def get_mock(dm_sim):
    """
    load dark matter simulation halo catalogue and make completeness cuts
    """
    if dm_sim=='Bolshoi_250':
        #open Bolshoi halo catalogue
        box_size=250
        filepath = cu.get_output_path() + 'processed_data/Multidark/Bolshoi/halo_catalogues/'
        halo_catalogue = 'hlist_1.00030.list'
        f =  h5py.File(filepath+halo_catalogue+'.hdf5', 'r')
        HC = f.get(halo_catalogue) #halo catalogue
        #for name in HC.dtype.names: print(name)
    elif dm_sim=='Chinchilla_125':
        #open Chinchilla halo catalogue
        box_size=125
        filepath = cu.get_output_path() + 'processed_data/Chinchilla/halo_catalogues/Lb125/'
        halo_catalogue = 'hlist_1.00000.list'
        f =  h5py.File(filepath+halo_catalogue+'.hdf5', 'r')
        HC = f.get(halo_catalogue) #halo catalogue
        #for name in HC.dtype.names: print(name)
    elif dm_sim=='Chinchilla_250':
        #open Chinchilla halo catalogue
        box_size=125
        filepath = cu.get_output_path() + 'processed_data/Chinchilla/halo_catalogues/Lb250/'
        halo_catalogue = 'hlist_1.00000.list'
        f =  h5py.File(filepath+halo_catalogue+'.hdf5', 'r')
        HC = f.get(halo_catalogue) #halo catalogue
        #for name in HC.dtype.names: print(name)
    else: print('simulation not available.')
    
    if dm_sim=='Bolshoi_250':
        mp_bolshoi = 1.35e8
        mpeak_cut = mp_bolshoi*200.0
        keep = (HC['Mpeak']>mpeak_cut)
        HC = HC[keep]
        #keep = (HC['vmax']>50.0)
        #HC = HC[keep]
        #print("number of (sub-)haloes in catalogue after cuts: {0}".format(len(HC)))
    if dm_sim=='Chinchilla_125':
        mp_bolshoi = 1.35e8
        mp_chinchilla = mp_bolshoi/8.0
        mpeak_cut = mp_chinchilla*200.0
        keep = (HC['Mpeak']>mpeak_cut)
        HC = HC[keep]
        #keep = (HC['vmax']>50.0)
        #HC = HC[keep]
        #print("number of (sub-)haloes in catalogue after cuts: {0}".format(len(HC)))
    if dm_sim=='Chinchilla_250':
        mp_bolshoi = 1.35e8
        mpeak_cut = mp_bolshoi*200.0
        keep = (HC['Mpeak']>mpeak_cut)
        HC = HC[keep]
        #keep = (HC['vmax']>50.0)
        #HC = HC[keep]
        #print("number of (sub-)haloes in catalogue after cuts: {0}".format(len(HC)))
     
    return HC, box_size


def get_sdss_sample():
    
    from scipy import interpolate
    
    filepath = cu.get_output_path() + 'processed_data/NYU_VAGC/'
    galaxy_catalogue = 'nyu_lss_mpa_vagc_dr7'
    f =  h5py.File(filepath+galaxy_catalogue+'.hdf5', 'r')
    GC = f.get(galaxy_catalogue) #halo catalogue
    GC = np.array(GC)
    #for name in GC.dtype.names: print(name)
    
    #trim the catalogue a bit
    zmin = 0.01
    zmax = 0.2
    selection = (GC['M']<=17.6) & (GC['Z']<zmax) & (GC['Z']>zmin) &\
                (GC['ZTYPE']==1) & (GC['FGOTMAIN']>0)
    GC = GC[selection]
    
    sm_key = 'sm_MEDIAN'
    GC[sm_key] = GC[sm_key]+np.log10(0.7**2.0)
    
    #make cuts on data to get a clean and complete sample
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    z = np.linspace(0.001,1,1000)
    dL = cosmo.luminosity_distance(z).value
    
    #make cheater fit to dl(z) function
    dL_z = interpolate.interp1d(z,dL,kind='linear')
    Mstar_lim = (4.852 + 2.246*np.log10(dL) + 1.123*np.log10(1+z) - 1.186*z)/(1.0-0.067*z)
    
    #dL = cosmo.luminosity_distance(GC['Z']).value
    dL = dL_z(GC['Z'])
    z = GC['Z']
    LHS = (4.852 + 2.246*np.log10(dL) + 1.123*np.log10(1.0+z) - 1.186*z)/(1.0-0.067*z)
    keep = (GC[sm_key]>LHS)
    GC = GC[keep]
    
    return GC

if __name__ == '__main__':
    main()