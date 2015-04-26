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
from Data_Exp import *
import time 
import matplotlib.cm as mplcm
import matplotlib.colors as colors

def sorptionX_aw(RH,C,Xm) :
    a=RH/100
    return Xm*(C*a/(1.+C*a)+a/(1.-a))
    
def simul_num(Mat,Cond):
        
    Xinit=Mat.Mass.X_init*np.ones_like(x)
    nbsavet=int(Cond.Time.final / Cond.Time.interval_savet)+1
    
    SaveX=np.empty((nbsavet,Mat.Size.number))
    SaveX[0,:]=Xinit
    
    Fo=Mat.Mass.diffusivity*Cond.Time.step/dx**2
    Vect_m=-Fo*np.ones(Mat.Size.number-1)
    Vect_p=-Fo*np.ones(Mat.Size.number-1)
    Vect=1.+2*Fo*np.ones(Mat.Size.number)

    if Cond.BC.physics=='fixed' :
        Vect[0]=1.
        Vect_p[0]=0.
    elif Cond.BC.physics=='flux' :
        Vect[0]=1.+2*Fo
        Vect_p[0]=-2*Fo

    Vect_m[-1]=-1.
    Vect[-1]=1.
       
    alpha,beta=Ant.Thomas_alpha_beta(Vect_m,Vect,Vect_p)
 
    if Cond.BC.physics=='fixed' :
        Xm=0.03398551582055446
        CL=10.906023000635287
        HR=Cond.Gaz.humidity(time)
        Xsurf=sorptionX_aw(HR,CL,Xm)*Mat.Mass.DM/Mat.Size.volume

    n=0    
    X=Xinit

    for p in xrange(len(time)-1) :
        rhs=X
        if Cond.BC.physics=='fixed' :
            rhs[0]=Xsurf[p+1]
        elif Cond.BC.physics=='flux' :
            rhs[0]=rhs[0]+2*Cond.BC.values((p+1)*Cond.Time.step)\
                *1e-3*Cond.Time.step/Mat.Size.area/dx
                        
        rhs[-1]=0
        X=Ant.Thomas_x(Ant.Thomas_y(beta,rhs),alpha,Vect_p)
        if ((p+1)%int(Cond.Time.interval_savet/Cond.Time.step))==0 :
            n=n+1
            print "--------------------------------------"
            print "Iteration %d - Temps %.5g" %(p,(p+1)*Cond.Time.step)
            SaveX[n,:]=X
    
    
    savetime=np.arange(0,nbsavet*Cond.Time.interval_savet,Cond.Time.interval_savet)
    
    return savetime,SaveX

def create_colors(n,p=2,colormap='jet') :
    cm = plt.get_cmap(colormap)
    cNorm  = colors.PowerNorm(p,vmin=0, vmax=n-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    L=[scalarMap.to_rgba(i) for i in range(n)]
    return L

def catch_final_time(val,Cond):
    return np.floor(val/Cond.Time.interval_savet)*Cond.Time.interval_savet

def continuation():
    answer=raw_input("Do you want to continue [y/n] : ")            
    while answer.upper()== "N" :                
        print ""    
        print 30*"="
        print ""   
        print "Type Ctrl-C to exit"
        print ""
        print 30*"="
        print ""
        time.sleep(1)
        
if __name__ == '__main__':
    plt.close('all')
    print 30*"="
    print ""
    print "Projet Ethicon"
    print ""    
    print 30*"="
    print ""    
    
    file_in='Ethicon-4x4.cfg'
    file_cond='Conditions-4x4.cfg'

    Conditions = ImportData(file_cond)
    Material = ImportData(file_in)
    
    Material.Size.volume = Material.Size.length * Material.Size.width * Material.Size.height
    Material.Size.area = Material.Size.length * Material.Size.width 

#------------------------
#    FLAGS
    sauvegarde = True
    compare_exp = False
    typeExp = "DisCharge"
#------------------------

#------------------------
#    Experimental    
    if compare_exp == True :
        if typeExp=="Charge" :
            file_exp='Result_Sorption_brut.csv'
#            data_exp=loaddata(file_exp,(0,1)) 
#            data_exp=loaddata(file_exp,(2,3))  # mass_flux2
            data_exp=loaddata(file_exp,(4,5))   # mass_flux3 
#            data_exp=loaddata(file_exp,(6,7)) 
            Material.Mass.X_init=0
            print "Charge Experiment"            

        elif typeExp=="DisCharge" :            
            file_exp='Result_DisCharge.csv'
            data_exp=loaddata(file_exp,(0,1)) # mass_flux_discharge
            Material.Mass.X_init=(data_exp[0,-1]*1e-3-Material.Mass.DM)/Material.Size.volume
            print "Discharge Experiment"

        print "You will compare to experiment : ",file_exp
        Conditions.Time.final = catch_final_time(data_exp[-1,0],Conditions)
        print "Boundary condition : ",Conditions.BC.physics
        print "Source of flux in : ",Conditions.BC.values_filename
        
    else :
        print "Simulation with humidity as input"
        print "Boundary condition : ",Conditions.BC.physics
        print "Source of humity vs time: ",Conditions.Gaz.humidity_filename
        Conditions.Time.final = catch_final_time(Conditions.Gaz.humidity_param[-1,0],Conditions)
        
    continuation()
#------------------------

#------------------------                
# Simulation            
    x,dx=np.linspace(0,Material.Size.height,Material.Size.number,retstep=True)
    time=np.arange(0,Conditions.Time.final+Conditions.Time.step,Conditions.Time.step)
    
    savetime,SaveX = simul_num(Material,Conditions)
#------------------------

#------------------------
# Plots    
    plt.figure(1)
    plt.xlabel(r"Distance $\left[mm\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3} \right]$")
    plt.grid(True)
    Lcolors=create_colors(len(savetime),p=0.1,colormap='jet')
    for i,c in zip(range(len(savetime)),Lcolors):
        plt.plot(x*1e3,SaveX.T[:,i],color=c)
    head_title=r"Water content profiles every "+str(int(Conditions.Time.interval_savet))+" s"
    head_title=head_title+"\n Mass diffusivity : D = "+str(Material.Mass.diffusivity)+r" $m^{2} \cdot s^{-1}$"
    plt.title(head_title)
    plt.grid(True)

    plt.figure(2)
    if compare_exp and typeExp=="Charge" :
        plt.plot(data_exp[:,0]/3600,(data_exp[:,1]-data_exp[0,1])*1e-3/Material.Size.volume,ls='None',\
            marker='o',mec='blue',mfc='white',label=r'$\bar{C}_{exp}$')
    elif compare_exp and typeExp=="DisCharge" :
        plt.plot(data_exp[:,0]/3600,(data_exp[:,1]*1e-3-Material.Mass.DM)/Material.Size.volume,ls='None',\
            marker='o',mec='blue',mfc='white',label=r'$\bar{C}_{exp}$')
        
    plt.plot(savetime/3600,SaveX.mean(axis=1),'k-',linewidth=3,label=r"$\bar{C}^{simul}$")
    plt.plot(savetime/3600,SaveX[:,0],'g-',label=r'$C^{simul}_{surf}$')
    plt.plot(savetime/3600,SaveX[:,-1],'r-',label=r'$C^{simul}_{deep}$')
    plt.xlabel(r"Time $\left[h\right]$")
    plt.ylabel(r"Water\ Content $\left[kg \cdot m^{-3}\right]$")
    plt.title(r"Water content vs. time")
    plt.legend()
    plt.grid(True)
#------------------------
    
#------------------------
# Sauvegarde en PDF        
    directory="./Resultats/"
    
    if sauvegarde==True :
        plt.figure(1)
        fname1="Figure_1_"+str(Material.Mass.diffusivity)+".pdf"
        plt.savefig(directory+fname1)
        plt.figure(2)
        fname2="Figure_2_"+str(Material.Mass.diffusivity)+".pdf"
        plt.savefig(directory+fname2)
#------------------------
        
