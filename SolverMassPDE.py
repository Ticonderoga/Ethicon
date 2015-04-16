#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
 Created on Mon Jan 19 15:02:44 2015

 Author : Philippe Baucour philippe.baucour@univ-fcomte.fr
"""
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,MA 02110-1301,USA.


import Routines.AnToolsPyx as Ant
import Routines.DiagHumid as Diag
from Routines.ImportData import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as scint

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
def sorptionX_aw(a,C,Xm) :
    return Xm*(C*a/(1.+C*a)+a/(1.-a))
    
def simul_num(Mat,Cond):
        
    Xinit=Mat.Mass.X_init*np.ones_like(x)
    nbsavet=int(Cond.Time.final / Cond.Time.interval_savet)+1
    
    SaveX=np.empty((nbsavet,Mat.Size.number))
    SaveX[0,:]=Xinit
    
    Fo=Mat.Mass.diffusivity*Cond.Time.step/dx**2
    Vect_m=-Fo*np.ones(Mat.Size.number-1)
    Vect_p=-Fo*np.ones(Mat.Size.number-1)
    Vect=1+2*Fo*np.ones(Mat.Size.number)
    Vect[0]=1.
    Vect_p[0]=0.
    Vect_m[-1]=-1.
    Vect[-1]=1.
    
    alpha,beta=Ant.Thomas_alpha_beta(Vect_m,Vect,Vect_p)
#    HR=Cond.Gaz.humidity(time)
    Xsurf=sorptionX_aw(HRsurf*np.ones_like(time),CL,Xm)
    n=0    
    X=Xinit
    for p in xrange(len(time)-1) :
        rhs=X
        rhs[0]=Xsurf[p+1]
        rhs[-1]=0
        X=Ant.Thomas_x(Ant.Thomas_y(beta,rhs),alpha,Vect_p)
        if ((p+1)%int(Cond.Time.interval_savet/Cond.Time.step))==0 :
            n=n+1
            print "======================"
            print "Iteration %d - Temps %.5g" %(p,(p+1)*Cond.Time.step)
            SaveX[n,:]=X
    
    
    savetime=np.arange(0,nbsavet*Cond.Time.interval_savet,Cond.Time.interval_savet)
    
    return savetime,SaveX

def loaddata(filename,columns) :
    ts2 = np.genfromtxt(filename, delimiter=",", filling_values=np.nan,skiprows=1,usecols=columns)
    ts2=ts2[np.isfinite(ts2)]
    ts2=np.reshape(ts2,(-1,2))
    return ts2

def interpdata(ts,timeint,derivation=0) :
    new_value_ts=scint.interp1d(ts[:,0],ts[:,1],bounds_error=False,fill_value=np.nan)(timeint)
    index=np.isfinite(new_value_ts)
    ts=np.c_[timeint[index],new_value_ts[index]]
    mass_flux=savitzky_golay(ts[:,1],window_size=21,order=3,deriv=derivation)/np.gradient(ts[:,0])
    return mass_flux

if __name__ == '__main__':
    
    file_in='Ethicon.cfg'
    file_cond='Conditions.cfg'
    file_exp='Result_Charge.csv'

    ts1=loaddata(file_exp,(0,1)) 
    ts2=loaddata(file_exp,(2,3)) 
    ts3=loaddata(file_exp,(4,5)) 
    ts4=loaddata(file_exp,(6,7)) 

    time_final_ts1=ts1[-1,0]
    time_ts1=np.arange(0,time_final_ts1,600.)
    mass_flux1=interpdata(ts1,time_ts1,derivation=1) 
    np.savetxt('mass_flux1.dat',np.c_[time_ts1,mass_flux1],header="time [s],mass_flux")
    
    time_final_ts2=ts2[-1,0]
    time_ts2=np.arange(0,time_final_ts2,600.)    
    mass_flux2=interpdata(ts2,time_ts2,derivation=1)
    np.savetxt('mass_flux2.dat',np.c_[time_ts2,mass_flux2],header="time [s],mass_flux")

    time_final_ts3=ts3[-1,0]
    time_ts3=np.arange(0,time_final_ts3,600.)
    mass_flux3=interpdata(ts3,time_ts3,derivation=1) 
    np.savetxt('mass_flux3.dat',np.c_[time_ts3,mass_flux3],header="time [s],mass_flux")

    time_final_ts4=ts4[-1,0]
    time_ts4=np.arange(0,time_final_ts4,600.)
    mass_flux4=interpdata(ts4,time_ts4,derivation=1) 
    np.savetxt('mass_flux4.dat',np.c_[time_ts4,mass_flux4],header="time [s],mass_flux")
    
    

    plt.close('all')
    print 30*"="
    print ""
    print "Projet Ethicon"
    print ""    
    print 30*"="
    plt.close('all')
    Material = ImportData(file_in)
    Conditions = ImportData(file_cond)
    
    x,dx=np.linspace(0,Material.Size.height,Material.Size.number,retstep=True)
    time=np.arange(0,Conditions.Time.final+Conditions.Time.step,Conditions.Time.step)



    
    
    Xm=0.03398551582055446
    CL=10.906023000635287
#    Xfinal = ts2.values[-1]
#    fopt= lambda x : Xfinal-sorptionX_aw(x,CL,Xm)
#    HRsurf=sp.optimize.fsolve(fopt,0.2)[0]
    savetime,SaveX = simul_num(Material,Conditions)





    plt.figure(1)
    plt.xlabel(r"Distance $\left[mm\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3} \right]$")
    plt.grid(True)
    plt.plot(x*1e3,SaveX.T)
#    plt.legend()
    plt.grid(True)
    plt.figure(2)
    plt.plot(savetime/3600,SaveX.mean(axis=1),label=r"X_{mean}")
    plt.plot(savetime/3600,SaveX[:,0],label=r'$X_{surf}$')
    plt.xlabel(r"Time $\left[h\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3}\ or\ \%\right]$")
    plt.legend()
    plt.grid(True)
#    plt.plot(ts2.index/3600,ts2.values)
