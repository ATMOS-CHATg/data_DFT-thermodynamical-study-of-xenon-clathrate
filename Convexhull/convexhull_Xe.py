#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:19:15 2024

@author: meko
"""
from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import ConvexHull
from scipy.optimize import curve_fit
import os
import re


plt.rcParams.update({
    'axes.linewidth': 2,
    'axes.grid': False,
    'grid.linestyle': '--',
    'font.size': 14,
})

script_dir = os.path.dirname(__file__)
print(script_dir)

xcfunc = "PBE"  # choice of  functional
guest = "Xe"         # choice of the guest molecule
sys = "sI"   # choice of the hydrate phase

def get_numbers_from_filename(filename):
    return re.findall(r'\d+', filename)[-1]


def birch_murnaghan_fit(a,E0,a0,B0,dB0):
    term1=9*a0**3*B0/16
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
    return (birch_murnaghan_fit(a,E0,a0,B0,dB0)+pressure(a,E0,a0,B0,dB0)*a**3)

p=1.0
phase=["XI","II","XV","VIII"]           # find the reference value for ice with pressure P
interpice = {}
listint=[]
for i in phase:
    dataice = np.loadtxt(f"{script_dir}/ICE/{xcfunc}/{i}/data-{xcfunc}-{i}.dat",usecols=[2,4])
    interpice[i] = interpolate.interp1d(dataice[:,0],dataice[:,1],kind='nearest',fill_value="extrapolate")
    listint.append(interpice[i](p))
refice = min(listint)

print(listint)
phaseguest=["FCC","HCP"]           # find the reference value for guest with pressure P
interpgas = {}
listgint=[]
for i in phaseguest:
    datagas = np.loadtxt(f"{script_dir}/{guest}/{xcfunc}/{i}/data-{xcfunc}-{i}.dat",usecols=[2,4])
    interpgas[i] = interpolate.interp1d(datagas[:,0],datagas[:,1],kind='nearest',fill_value="extrapolate")
    listgint.append(interpgas[i](p))
refgas = min(listgint)


pointsTOT = [(0,0),(1,0)]
fig, ax = plt.subplots(figsize=(10,3.5))
points = []   
nwater = 46
    

for filename in os.listdir(f"{script_dir}/Clathrate_{guest}/{xcfunc}/"):
    print(filename)
    ngas = int(get_numbers_from_filename(filename))  # get the composition from file name
    ntot = ngas+nwater
    comp = nwater/ntot   
    
    # interpolation of the clathrate data        
                                    
    dataclath=np.loadtxt(f"{script_dir}/Clathrate_{guest}/{xcfunc}/"+filename) 
    adft = dataclath[:,0]
    edft = dataclath[:,1]
    popt,pcov = curve_fit(birch_murnaghan_fit,adft,edft,p0=[min(edft),np.mean(adft),5,0.1])
    xpre = np.linspace(min(adft)-3,popt[1],1000)
    datasave=np.array([pressure(xpre,*popt)*160.2,enthalpy(xpre,*popt)])
    refint = interpolate.interp1d(datasave[0],datasave[1],kind='nearest',fill_value="extrapolate")
    print(refint(p))
    
    # relative enthalpy calculation at a pressure p      
                                            
    enthrel = (refint(p)-nwater*refice-ngas*refgas)/ntot   # non bonding energy
    print(enthrel)
    points.append([ngas,comp,enthrel])
    pointsTOT.append([comp,enthrel])
points = np.array(points)
print(points[:,1],points[:,2])


plt.plot(points[:,1],points[:,2],'o') 
pointsTOT = np.array(pointsTOT)

# Place differnt points of relative enthalpy pon a convexhull

hull = ConvexHull(pointsTOT)
for simplex in hull.simplices[0:]:
    if pointsTOT[simplex[0],1] < 0:
        plt.plot(pointsTOT[simplex,0],pointsTOT[simplex,1],'r--')
    plt.title("P="+str(p)+" GPa", fontname='Garamond',fontsize=18, fontweight='bold')
    plt.xlim(0,1)
    plt.ylim(-0.03, 0.05)
    plt.xlabel('$\mathregular{Fraction\ {[H_{2}O]}/([{H_{2}O] + [Xe]})}$', fontsize=18, fontname='Garamond',fontweight='bold')
    plt.ylabel("$\Delta$H (eV/molecule)", fontsize=18, fontname='Garamond',fontweight='bold')
    plt.xticks(fontsize=18, fontname='Garamond',fontweight='bold')
    plt.yticks(fontsize=18, fontname='Garamond',fontweight='bold')
plt.subplots_adjust(left=0.300, right=0.780, bottom=0.185, top=0.885, wspace=0.2, hspace=0.2)
plt.savefig(f"{script_dir}/Data-{xcfunc}_{p}GPa.pdf",bbox_inches='tight', pad_inches=0.05)
plt.show()


