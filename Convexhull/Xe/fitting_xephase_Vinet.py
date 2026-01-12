#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:12:31 2024

@author: meko
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

script_dir = os.path.dirname(__file__)

#################################################################################
#Definition of the phase under study, the volume coefficient, and the number of atoms

xcfunc = "PBE"        
sys = "HCP"
if sys =="FCC":
   coeffvol=1
   nbxenon=4
elif sys =="HCP":
   coeffvol=1.414
   nbxenon=2
##################################################################################"

# Define the Vinet equation
def Vinet_fit(a,E0,a0,B0,dB0):
    term1=B0*a0**3*coeffvol/(dB0-1)**2
    term2=(dB0-1)
    term3=((a/a0)-1)
    return E0 + 4*term1-2*term1*(3*term2*term3+2)*np.exp(-1.5*term2*term3)

def pressure(a,E0,a0,B0,dB0):
     term1=3*B0*((a/a0))**(-2)
     term2=(1-(a/a0))
     term3=(dB0-1)
     return term1*term2*np.exp(1.5*term3*term2)

def enthalpy(a,E0,a0,B0,dB0):
    return (Vinet_fit(a,E0,a0,B0,dB0)+pressure(a,E0,a0,B0,dB0)*a**3*coeffvol)/nbxenon
###############################################################################
# Load data from file

data = np.loadtxt(f"{script_dir}/{xcfunc}/{sys}/energydft.dat")
a= data[:,0]
E_values= data[:,1]



popt,pcov = curve_fit(Vinet_fit,a,E_values,p0=[min(E_values),np.mean(a),1,.1], maxfev=1000)
perr=np.sqrt(np.diag(pcov))
xpre = np.linspace(min(a)-0.25,popt[1],10000)

###############################################################################
# Saving data as : a,V,P,E,H

datasave=np.array([xpre,xpre**3*coeffvol/nbxenon,pressure(xpre,*popt)*160.2,Vinet_fit(xpre,*popt)/nbxenon,enthalpy(xpre,*popt)])
dataent=[xpre,enthalpy(xpre,*popt)]
np.savetxt(f"{script_dir}/{xcfunc}/{sys}/data-{xcfunc}-{sys}.dat",np.column_stack(datasave))

################################################################################
plt.gcf().subplots_adjust(wspace=0.5)
plt.subplot(221)
plt.plot(a,Vinet_fit(a,popt[0],popt[1], popt[2],popt[3]),color='r',label='fitting curve' )
plt.plot(a,E_values,'o',color='g', label='DFT') 
plt.ylabel('Potential energy(eV)')
plt.title (f"fit Xe_{xcfunc}")
plt.subplot(223)
plt.xlabel('Lattice parameter (Ang) ')
plt.ylabel('Pressure(GPa)')
plt.plot(xpre,pressure(xpre,*popt)*160.2,color='r')

plt.subplot(222)
plt.ylabel('Potential energy (eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,Vinet_fit(xpre,*popt)/nbxenon, color='r')
plt.subplot(224)
plt.xlabel('Pression(GPa) ')
plt.ylabel('Enthalpy(eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,enthalpy(xpre,*popt), color='r')
plt.legend()
plt.savefig(f"{script_dir}/{xcfunc}/{sys}/properties-{xcfunc}-{sys}.pdf")
plt.show()
print(popt)
print("\nE0=",popt[0],"±",perr[0],"eV")
print("a0=",popt[1],"±",perr[1],"Ang")
print("B0 =", popt[2]*160.2,"±",perr[2]*160.2,"GPa")  
print("dB0=",popt[3],"±",perr[3])

    
    
    
    
    
    
