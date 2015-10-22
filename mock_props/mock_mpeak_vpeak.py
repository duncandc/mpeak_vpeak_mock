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
from mpl_toolkits.mplot3d import Axes3D

def main():

    #open mock
    catalogue = 'sm_9.5_s0.1_sfr_c-1.0_250'
    filepath = cu.get_output_path() + 'processed_data/campbell_mocks/'
    f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f.get(catalogue)
    mock = np.array(mock)
    print(mock.dtype.names)
    
    host = (mock['upid']==-1)
    sub = (mock['upid']!=-1)
    
    LHS = -11.0
    blue = (mock['SSFR']>LHS)
    red = (mock['SSFR']<LHS)
    
    bins = np.arange(0,3,0.01)
    bins = 10.0**bins
    bin_centers = (bins[:-1]+bins[1:])/2.0
    bin_widths = bins[1:]-bins[:-1]
    left,right = bins[:-1],bins[1:]
    X = np.array([left,right]).T.flatten()
    selection = (mock['Mpeak']>10**12) & (mock['Mpeak']<10**12.5)
    counts_cen = np.histogram(mock['Vpeak'][selection&host],bins=bins)[0]
    counts_sat = np.histogram(mock['Vpeak'][selection&sub],bins=bins)[0]
    counts_cen = counts_cen/(np.sum(counts_cen)*1.0)
    counts_sat = counts_sat/(np.sum(counts_sat)*1.0)
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    p1 = plt.bar(bins[:-1], counts_cen, width=bin_widths, color='green', alpha=0.5, lw=0)
    Y = np.array([counts_cen,counts_cen]).T.flatten()
    plt.plot(X,Y,color='green')
    p2 = plt.bar(bins[:-1], counts_sat, width=bin_widths, color='orange', alpha=0.5, lw=0)
    Y = np.array([counts_sat,counts_sat]).T.flatten()
    plt.plot(X,Y,color='orange')
    plt.xlim([100,400])
    plt.xlabel(r'$V_{\rm peak} [{\rm km/s}]$')
    plt.ylabel(r'$P(V_{\rm peak})$')
    plt.title(r'$12.0<\log(M_{\rm peak}/h^{-1}M_{\odot})<12.5$')
    plt.legend((p1,p2),('host','sub'),fontsize=10.0,frameon=False)
    plt.show()
    
    bins = np.arange(0,3,0.01)
    bins = 10.0**bins
    bin_centers = (bins[:-1]+bins[1:])/2.0
    bin_widths = bins[1:]-bins[:-1]
    left,right = bins[:-1],bins[1:]
    X = np.array([left,right]).T.flatten()
    selection = (mock['Mpeak']>10**13) & (mock['Mpeak']<10**13.5)
    counts_cen = np.histogram(mock['Vpeak'][selection&host],bins=bins)[0]
    counts_sat = np.histogram(mock['Vpeak'][selection&sub],bins=bins)[0]
    counts_cen = counts_cen/(np.sum(counts_cen)*1.0)
    counts_sat = counts_sat/(np.sum(counts_sat)*1.0)
    counts_cen = counts_cen/(np.sum(counts_cen))
    counts_sat = counts_sat/(np.sum(counts_sat))

    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    plt.bar(bins[:-1], counts_cen, width=bin_widths, color='green', alpha=0.5, lw=0)
    Y = np.array([counts_cen,counts_cen]).T.flatten()
    plt.plot(X,Y,color='green')
    plt.bar(bins[:-1], counts_sat, width=bin_widths, color='orange', alpha=0.5, lw=0)
    Y = np.array([counts_sat,counts_sat]).T.flatten()
    plt.plot(X,Y,color='orange')
    plt.xlim([200,700])
    plt.xlabel(r'$V_{\rm peak}~[{\rm km/s}]$')
    plt.ylabel(r'$P(V_{\rm peak})$')
    plt.show()
    
    
if __name__ == '__main__':
    main()