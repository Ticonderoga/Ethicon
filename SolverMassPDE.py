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
import scipy.interpolate as scinterp
    
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
    Xsurf=sorptionX_aw(0.2*np.ones_like(time),CL,Xm)
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


if __name__ == '__main__':
    
    file_in='Ethicon.cfg'
    file_cond='Conditions.cfg'

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
