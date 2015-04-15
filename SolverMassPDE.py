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
    
    Cinit=Material.Mass.Conc_init*np.ones_like(x)
    nbsavet=int(Conditions.Time.final / Conditions.Time.interval_savet)+1
    
    SaveC=np.empty((nbsavet,Material.Size.number))
    SaveC[0,:]=Cinit
    
    Fo=Material.Mass.diffusivity*Conditions.Time.step/dx**2
    Vect_m=-Fo*np.ones(Material.Size.number-1)
    Vect_p=-Fo*np.ones(Material.Size.number-1)
    Vect=1+2*Fo*np.ones(Material.Size.number)
    Vect[0]=1.
    Vect_p[0]=0.
    Vect_m[-1]=-1.
    Vect[-1]=1.
    
    alpha,beta=Ant.Thomas_alpha_beta(Vect_m,Vect,Vect_p)
    Csurf=Conditions.Gaz.humidity(time)
    n=0    
    C=Cinit
    for p in xrange(len(time)-1) :
        rhs=C
        rhs[0]=Csurf[p+1]
        rhs[-1]=0
        C=Ant.Thomas_x(Ant.Thomas_y(beta,rhs),alpha,Vect_p)
        if ((p+1)%int(Conditions.Time.interval_savet/Conditions.Time.step))==0 :
            n=n+1
            print "======================"
            print "Iteration %d - Temps %.5g" %(p,(p+1)*Conditions.Time.step)
            SaveC[n,:]=C
    
    
    savetime=np.arange(0,nbsavet*Conditions.Time.interval_savet,Conditions.Time.interval_savet)
    plt.figure(1)    
    plt.xlabel(r"Distance $\left[mm\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3} \right]$")
    plt.grid(True)
    plt.plot(x*1e3,SaveC.T)
#    plt.legend()
    plt.grid(True)
    plt.figure(2)
    plt.plot(savetime/3600,SaveC.mean(axis=1),label=r"C_{mean}")
    plt.plot(time/3600,Csurf,label=r'RH')
    plt.xlabel(r"Time $\left[h\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3}\ or\ \%\right]$")
    plt.legend()
    plt.grid(True)