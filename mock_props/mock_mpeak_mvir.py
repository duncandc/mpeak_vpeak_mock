#Duncan Campbell
#March 2015
#Yale University
#examine the relation between mpeak and mvir

#load packages
from __future__ import print_function, division
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    rho = -1.0
    sigma = 0.2

    #open halo catalogue (used to account for empty haloes)
    sys.path.insert(0, '../mocks/')
    from make_mock import get_mock
    mock, Lbox = get_mock('Bolshoi')
    mock = np.array(mock)
    print(mock.dtype.names)
    
    hosts = (mock['upid']==-1)
    subs = (mock['upid']!=-1)
    
    never_sub = ((mock['First_Acc_Scale']==1.0) & (hosts))
    
    print(np.sum(hosts),np.sum(never_sub))

    sys.exit()

    fig = plt.figure()
    plt.plot(mock['Mpeak'][hosts], mock['mvir'][hosts], '.', ms=2, alpha=0.1, color='green')
    plt.plot(mock['Mpeak'][never_sub], mock['mvir'][never_sub], '.', ms=2, alpha=0.1, color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.show(block=True)
    
    
    fig = plt.figure()
    plt.plot(mock['Mpeak'][subs], mock['mvir'][subs], '.', ms=2, alpha=0.1)
    plt.xscale('log')
    plt.yscale('log')
    plt.show(block=False)
    
    bins=np.arange(10,15,0.1)
    inds = np.digitize(np.log10(mock['Mpeak']),bins=bins)
    
    sub_inds = (inds==15)
    host_counts, bins = np.histogram(np.log10(mock['mvir'][sub_inds & hosts]),bins)
    sub_counts, bins = np.histogram(np.log10(mock['mvir'][sub_inds & subs]),bins)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    bin_widths = bins[1:]-bins[:-1]
    left,right = bins[:-1],bins[1:]
    X = np.array([left,right]).T.flatten()
    
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2, top=0.9)
    p1 = plt.bar(bins[:-1], host_counts, width=bin_widths, color='green', alpha=0.5, lw=0)
    Y = np.array([host_counts,host_counts]).T.flatten()
    p1 = plt.bar(bins[:-1], sub_counts, width=bin_widths, color='orange', alpha=0.5, lw=0)
    Y = np.array([sub_counts,sub_counts]).T.flatten()
    plt.plot(X,Y,color='orange')
    plt.show()
    

if __name__ == '__main__':
    main()