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
xcfunc = "revPBE"  

#################################################################################
#Definition of the phase under study, the volume coefficient, and the number of atoms

sys = "FCC"
if sys =="FCC":
   coeffvol=1
   nbxenon=4
elif sys =="HCP":
   coeffvol=1.414
   nbxenon=2
##################################################################################"
# Define the Birch-Murnaghan equation
def birch_murnaghan_fit(a,E0,a0,B0,dB0):
    term1=9*a0**3*coeffvol*B0/16
    term2=((a0/a)**2-1)**3
    term3=((a0/a)**2-1)**2
    term4=(6-4*(a0/a)**2)
    return E0 + term1*(term2*dB0 + term3*term4)

def pressure(a,E0,a0,B0,dB0):
    term1=((a0/a)**7-(a0/a)**5)
    term2=(3*(dB0-4)/4)
    term3=((a0/a)**2-1)
    return 3*B0/2*term1*(1+term2*(term3))
def enthalpy(a,E0,a0,B0,dB0):
    return (birch_murnaghan_fit(a,E0,a0,B0,dB0)+pressure(a,E0,a0,B0,dB0)*a**3*coeffvol)/nbxenon
# Load data from file
data = np.loadtxt(f"{script_dir}/{xcfunc}/{sys}/energydft.dat")
a= data[::,0]
E_values= data[::,1]


mean_a=sum(a)/len(a)
popt,pcov = curve_fit(birch_murnaghan_fit,a,E_values,p0=[min(E_values),mean_a,1,0.1])
perr=np.sqrt(np.diag(pcov)) 
xpre = np.linspace(min(a),popt[1],10000)

###############################################################################
# Saving data as : a,V,P,E,H

datasave=np.array([xpre,xpre**3*coeffvol/nbxenon,pressure(xpre,*popt)*160.2,birch_murnaghan_fit(xpre,*popt)/nbxenon,enthalpy(xpre,*popt)])
dataent=[xpre,enthalpy(xpre,*popt)]
(f"{script_dir}/{xcfunc}/{sys}/data-{xcfunc}-{sys}.dat",np.column_stack(datasave))

################################################################################
plt.gcf().subplots_adjust(wspace=0.5)
plt.subplot(221)
plt.plot(a,birch_murnaghan_fit(a,popt[0],popt[1], popt[2],popt[3]),color='r',label='fitting curve' )
plt.plot(a,E_values,'o',color='g', label='DFT') 
plt.ylabel('Potential energy(eV)')
plt.title ('fit Xe_vdWDF2')


plt.subplot(223)
plt.xlabel('Lattice parameter (Ang) ')
plt.ylabel('Pressure(GPa)')
plt.plot(xpre,pressure(xpre,*popt)*160.2,color='r')

plt.subplot(222)
plt.ylabel('Potential energy (eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,birch_murnaghan_fit(xpre,*popt)/nbxenon, color='r')
plt.subplot(224)
plt.xlabel('Pression(GPa) ')
plt.ylabel('Enthalpy(eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,enthalpy(xpre,*popt), color='r')
plt.legend()
plt.show()
plt.savefig(f"{script_dir}/{xcfunc}/{sys}/properties-{xcfunc}-{sys}.pdf")


print("\nE0=",popt[0],"±",perr[0],"eV")
print("a0=",popt[1],"±",perr[1],"Ang")
print("B0 =", popt[2],"±",perr[2],"GPa") 
print("dB0=",popt[3],"±",perr[3])

    
    
    
    
    
    
