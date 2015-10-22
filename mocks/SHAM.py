#Duncan Campbell
#Yale University
#March 2015

"""
Build a SHAM mock.
"""

#load modules
from __future__ import print_function, division
import numpy as np
from scipy import interpolate
from scipy import optimize
from scipy import integrate
from scipy.misc import derivative
import matplotlib.pyplot as plt
import sys

__all__=['AM', 'make_SHAM_mock',\
         'abundance_function_from_tabulated', 'abundance_function_from_observations']


def AM(dn_dx, x, dn_dy, y, P_x, tol=0.0001):
    """
    Abundance match galaxy property 'x' to halo property 'y'.
    
    Functional forms of the abundance of galaxies and haloes are required.
    
    Determines P(x_gal | y_halo), given the centered distribution, i.e. all other moments.
    
    In detail P should be of the form: P_x(y, mu_xy(y)), where y is a property of the halo
    and mu_xy(y) is the mean of the distribution as a function of y.  This function simply
    solves for the form of mu_xy, returning the full conditional probability function.
    
    Parameters
    ==========
    dn_dx: function
        galaxy abundance as a function of property 'x'
    
    dn_dy: function
        halo abundance as a function of property 'y'
    
    P_x: function
        The central moment of P(x_gal | y_halo)
    
    xbins: array_like
        bins to sample the galaxy abundance function.  These are used during numerical 
        integrations, so the sampling should appropriately sample the galaxy abundance 
        function 
    
    ybins: array_like
        bins to sample the halo abundance function  These are used during numerical 
        integrations, so the sampling should appropriately sample the halo abundance 
        function 
    
    tol: float, optional
        tolerance of error on of the 1st moment of P(x| y).
    
    Returns
    =======
    P(x| y): function 
    """
    
    x = np.array(x)
    y = np.array(y)
    dn_dx = np.array(dn_dx)
    dn_dy = np.array(dn_dy)
    
    #check numerical values of dn_dx and dn_dy
    #trim down ranges if they return numbers that are too small
    keep = (dn_dx>10.0**(-20.0))
    if not np.all(keep):
        print("Triming x-range to keep number densities above 1e-20")
        x = x[keep]
        dn_dx = dn_dx[keep]
    keep = (dn_dy>10.0**(-20.0))
    if not np.all(keep):
        print("Triming y-range to keep number densities above 1e-20")
        y = y[keep]
        dn_dy = dn_dy[keep]
    
    #convert tabulated abundance functions into function objects
    #galaxy abundance function
    dn_dx = _abundance_function_from_tabulated(x, dn_dx)
    #halo abundance function
    dn_dy = _abundance_function_from_tabulated(y, dn_dy)
    
    """
    #plot abundance functions
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(dn_dy.x,dn_dy(dn_dy.x))
    plt.plot(dn_dx.x,dn_dx(dn_dx.x))
    plt.yscale('log')
    plt.ylim([10**(-8),10])
    plt.xlim([8,16])
    plt.show()
    """
    
    #make sure bins are in sorted order
    x = np.sort(x)
    y = np.sort(y)
    
    #check for monotonicity
    #galaxy abundance function
    greater = np.array([False]*(len(x)-1))
    for i in range(0,len(x)-1):
        greater[i] = dn_dx(x[i])<dn_dx(x[i+1])
    if not (np.all(greater) | np.all(greater==False)):
        print("warning galaxy abundance function is not monotonic")
    if np.all(greater==False):
        x = x[::-1]
        inverse_x=True
    else: inverse_x=False
    #halo abundance function
    greater = np.array([False]*(len(y)-1))
    for i in range(0,len(y)-1):
        greater[i] = dn_dy(y[i])<dn_dy(y[i+1])
    if not (np.all(greater) | np.all(greater==False)):
        print("warning halo abundance function is not monotonic")
    Ngreater = np.sum(greater)
    Nless = len(greater) - Ngreater
    if Nless > Ngreater:
    #if np.all(greater==False):
        y = y[::-1]
        inverse_y=True
    else: inverse_y=False
    
    print("inverse_x",inverse_x)
    print("inverse_y",inverse_y)
    
    
    #calculate the cumulative abundance functions and their inverses
    #for haloes
    N_cum_halo = _cumulative_abundance(dn_dy, y, inverse=inverse_y)
    N_cum_halo_inv = interpolate.interp1d(N_cum_halo(y), y, kind='linear', bounds_error=False, fill_value=0.0)
    #for galaxies
    N_cum_gal = _cumulative_abundance(dn_dx, x, inverse=inverse_x)
    N_cum_gal_inv = interpolate.interp1d(N_cum_gal(x), x, kind='linear', bounds_error=False, fill_value=0.0)
    
    """
    #plot the cumulative abundances of galaxies and haloes
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(x,np.log10(N_cum_gal(x)),'.')
    plt.plot(y,np.log10(N_cum_halo(y)),'.')
    plt.xlabel(r'$x,y$')
    plt.ylabel(r'$N(>x,y)$')
    plt.xlim([8,16])
    plt.ylim([-8,1])
    plt.show(block=True)
    """
    
    ##########################################################
    #step 1: solve for first moment in the absence of scatter
    
    x_y1 = N_cum_gal_inv(N_cum_halo(y))
    x_y1 = interpolate.interp1d(y, x_y1, kind='linear', bounds_error=False, fill_value=0.0)
    
    """
    #plot x(y)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(y,x_y1(y),'.')
    plt.xlabel(r'$M_{\rm vir}$')
    plt.ylabel(r'$M_{*}$')
    plt.xlim([10,16])
    plt.ylim([9,12])
    plt.show(block=True)
    """
    
    #apply first estimate of the 1st moment
    P1 = lambda y: P_x(y, mu_xy=x_y1)
    
    ##########################################################
    #step 2: determine dNhalo/dx= integral[ dy * dNhalo/dy * (P(x|y) ]
    
    #define integrand
    def integrand(y,x):
        return P1(y).pdf(x)*dn_dy(y)
    
    """    
    #take a look at the form of the integral
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(y,integrand(y,8))
    plt.plot(y,integrand(y,9))
    plt.plot(y,integrand(y,10))
    plt.plot(y,integrand(y,11))
    plt.plot(y,integrand(y,11.5))
    plt.show()
    """
    
    dn_halo_dx = np.zeros(len(x)) #define array to store numeric integration result
    for i, xx in enumerate(x):
        f = lambda y: integrand(y,xx) #simplify the integrand
        dn_halo_dx[i] = integrate.simps(f(y),y)*-1.0
    
    #interpolate result to get a functional form of the result
    dn_halo_dx = interpolate.interp1d(x,dn_halo_dx,kind='linear', bounds_error=False, fill_value=0.0)
    #integrate to get cumulative distribution
    N_cum_halo_x = _cumulative_abundance(dn_halo_dx, x, inverse=inverse_x)
    N_cum_halo_x_inv = interpolate.interp1d(N_cum_halo_x(x), x, bounds_error=False, fill_value=0.0)
    
    """
    #plot N_cum_halo_x, N_cum_gal
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(x,N_cum_gal(x),color='blue')
    plt.plot(x,N_cum_halo_x(x),color='red')
    plt.xlabel(r'$M_{\rm *}$')
    plt.ylabel(r'$N(<M_{*})$')
    plt.yscale('log')
    plt.ylim([10**(-15),1])
    plt.xlim([10,16])
    plt.show(block=True)
    """
    
    #calculate x'(x)
    x2_x1 = N_cum_gal_inv(N_cum_halo_x(x))
    x2_x1 = interpolate.interp1d(x, x2_x1, kind='linear', bounds_error=False, fill_value=0.0)
    
    """
    #plot x(x) relationship
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(x,x2_x1(x))
    plt.show()
    """
    
    #solve for new mean
    x_y2 = x2_x1(x_y1(y))
    x_y2 = interpolate.interp1d(y, x_y2, kind='linear', bounds_error=False, fill_value=0.0)
    
    """
    #plot 1st two iterations of x(y)
    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(y,x_y1(y),color='red')
    plt.plot(y,x_y2(y),color='blue')
    plt.xlabel(r'$M_{\rm vir}$')
    plt.ylabel(r'$M_{*}$')
    plt.xlim([10,16])
    plt.ylim([8,12])
    plt.show(block=True)
    """
    
    #apply mean
    P2 = lambda y: P_x(y, mu_xy=x_y2)
    
    #check for convergence
    keep = (x_y2(y)>0.0)
    delta_xx = np.max(np.abs(x_y2(y)-x_y1(y))[keep]/x_y2(y)[keep])
    print("delta_xx:", delta_xx)
    if delta_xx<tol:
        return P2
    
    ##########################################################
    #step 3: iterate step 2 until convergence
    convergence=False
    x_yi = x_y2
    Pi = P2

    fig = plt.figure(figsize=(3.3,3.3))
    fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.9)
    plt.plot(y,x_y1(y),color='black')
    plt.plot(y,x_y2(y),color='grey')
    plt.xlabel(r'$M_{\rm vir}$')
    plt.ylabel(r'$M_{*}$')
    plt.xlim([10,16])
    plt.ylim([8,12])
    plt.show(block=False)
    
    while convergence==False:
        
        #define integrand
        def integrand(y,x):
            return Pi(y).pdf(x)*dn_dy(y)
    
        dn_halo_dx = np.zeros(len(x)) #define array to store numeric integration result
        for i, xx in enumerate(x):
            f = lambda y: integrand(y,xx) #simplify the integrand
            dn_halo_dx[i] = integrate.simps(f(y),y)*-1.0
    
        #interpolate result to get a functional form of the result
        dn_halo_dx = interpolate.interp1d(x,dn_halo_dx,kind='linear', bounds_error=False, fill_value=0.0)
        #integrate to get cumulative distribution
        N_cum_halo_x = _cumulative_abundance(dn_halo_dx, x, inverse=inverse_x)
        N_cum_halo_x_inv = interpolate.interp1d(N_cum_halo_x(x), x, bounds_error=False, fill_value=0.0)
    
        #calculate x'(x)
        x2_x1 = N_cum_gal_inv(N_cum_halo_x(x))
        x2_x1 = interpolate.interp1d(x, x2_x1, kind='linear', bounds_error=False, fill_value=0.0)
    
        #solve for new mean
        x_yii = x2_x1(x_yi(y))
        x_yii = interpolate.interp1d(y, x_yii, kind='linear', bounds_error=False, fill_value=0.0)
        
        Pi = lambda y: P_x(y, mu_xy=x_yii)
        
        keep = (x_yii(y)>0.0)
        
        delta_xx = np.max(np.abs(x_yii(y)[keep]-x_yi(y)[keep])/x_yii(y)[keep])
        ind = np.where((np.abs(x_yii(y)[keep]-x_yi(y)[keep])/x_yii(y)[keep])==delta_xx)[0]
        print("delta_xx:", delta_xx, y[keep][ind])
        
        #plot next two iterations
        plt.plot(y,x_yii(y),color='grey')
        plt.draw()
        
        if delta_xx<tol:
            convergence=True
        else:
            x_yi = x_yii
    
    plt.show(block=True)
    
    return Pi


def _cumulative_abundance(f,bins,inverse=True):
    """
    calculate the cumulative abundance function by integrating the differential abundance 
    function
    """
    
    from scipy import integrate
    from scipy import optimize
    
    N_sample = len(bins) #number of points to sample along function
    
    N0=10.0**(-20.0) #use a small number for the first value
    N0 = float(f(bins[0]))*np.fabs(bins[0]-bins[1])
    N_cum = integrate.cumtrapz(f(bins), bins, initial=N0)
    
    #if the integration range is reversed, the integral will be negative.  So fix this if 
    #this is the case.
    if inverse==True:
         N_cum = -1.0* N_cum
         N_cum[0] = -1.0* N_cum[0] #we set this by hand, so it was fine.
    
    #interpolate result to get a functional form of the result
    N_cum_func = interpolate.interp1d(bins, N_cum, kind='linear', bounds_error=False, fill_value=0.0)
    
    return N_cum_func


def _minimize_cum_dist_diff(fx, fy, xbins, ybins, tol=0.0001):
    """
    find the value which minimizes the difference between two function given a sampling of
    each function: xbins, ybins.
    """
    
    from scipy import optimize
    
    #define function to minimize
    def F(x, y):
        return np.fabs(fy(y)-fx(x))

    x_y = np.zeros(len(ybins)) #empty array to store result of minimization
    
    x0 = np.median(xbins) #initial guess for x in minimization
    
    bounds = [np.min(xbins),np.max(xbins)] #bounds to search for minimization over
    
    #keep the y value fixed and search for the x value
    for i,y in enumerate(ybins):
        f = lambda x: F(x, y) #simplify minimization function
        x_y[i] = optimize.minimize_scalar(f, bounds = bounds, method='bounded').x
    
    y_x = interpolate.interp1d(x_y, ybins, bounds_error=False, fill_value=0.0)
    x_y = interpolate.interp1d(ybins, x_y, bounds_error=False, fill_value=0.0)
    
    return x_y, y_x


def make_SHAM_mock(mock, P_xy, mock_prop='mvir', gal_prop='mstar'):
    """
    make a SHAM mock given a halo catalogue, mock, the probability of galaxy property 'x' 
    given a halo with property 'y', where mock[mock_prop] returns halo property 'y'.
    
    Parameters
    ==========
    mock: array_like
        structured array containing halo catalogue
    
    P_xy: function
        Function that returns x_gal given y_halo
    
    mock_prop: string
        key into mock which returns the halo property to build the SHAM mock
    
    Returns
    =======
    mock: structured array
        mock with new column containing galaxy property gal_prop
    """
    
    from numpy.lib.recfunctions import append_fields
    mock.appendfield(gal_prop)
    
    mock[gal_prop]=P_xy(mock[mock_prop])
    
    return mock


def _abundance_function_from_tabulated(x, dndx):
    """
    given a tabulated abundance function for property 'x', build a functional form.
    
    Parameters
    ==========
    x: array_like
        values of property 'x'
    
    dndx: array_like
        number density
    
    Returns
    =======
    dndx: function
    """
    from scipy import interpolate
    
    n_x = interpolate.interp1d(x, dndx, kind='linear')
    min_x = np.min(x) #minimum of the x range
    max_x = np.max(x) #maximum of the x range
    min_n = np.min(dndx) #minimum of the number density
    max_n = np.max(dndx) #maximum of the number density
    
    def f(x):
        return n_x(x)
    
    #store x and y ranges of the function
    setattr(f,'min_x',min_x)
    setattr(f,'max_x',max_x)
    setattr(f,'min_dn',min_n)
    setattr(f,'max_dn',max_n)
    
    #store tabulated form of the function
    setattr(f,'x',x)
    setattr(f,'dn',dndx)
    
    return f


def abundance_function_from_observations(x, weights, bins):
    """
    given observations of objects with property 'x' and weights, return an abundance 
    function
    
    Parameters
    ==========
    x: array_like
    
    weights: array_like
    
    bins: array_like
    
    Returns
    =======
    dndx: dn, x
    """
    
    if np.shape(weights)==():
        weights = np.array([weights]*len(x))
    
    n = np.histogram(x,bins,weights=weights)[0]
    bin_centers = (bins[:-1]+bins[1:])/2.0
    dx = bins[1:]-bins[:-1]
    dn = n/dx
    
    return dn, bin_centers
    

def _extrap1d(interpolator):
    
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def f(x):
        if np.shape(x)==():
            x = np.array([x])
        else: x = np.array(x)
        return np.array(map(pointwise, x))

    return f
    
    
