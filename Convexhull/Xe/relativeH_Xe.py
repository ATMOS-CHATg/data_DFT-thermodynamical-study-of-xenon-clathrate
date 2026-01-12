#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 11:46:46 2024

@author: meko
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os


plt.style.use('default')
plt.rcParams.update({
    'figure.figsize': (6, 6),
    'font.family': 'Garamond',
    'font.size': 14,
    'axes.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 0.8,
    'lines.markersize': 6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    "axes.labelsize": 20,
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True
})

script_dir = os.path.dirname(__file__)

xcfunc = "PBE"  # functional choice

def data(phase,ref):
    a = np.loadtxt(script_dir+"/"+xcfunc+"/"+phase+"/data-"+xcfunc+"-"+phase+".dat")

    return a[:,2],a[:,4]-ref(a[:,2])


sys=["FCC","HCP"]  

###############################################################################
# Calculation of relative enthalpies
dataref = np.loadtxt(f"{script_dir}/{xcfunc}/FCC/data-{xcfunc}-FCC.dat")
dataref = dataref[:,2],dataref[:,4]
datasave = dataref[0],dataref[1]-dataref[1]
refint = interpolate.interp1d(dataref[0],dataref[1],kind='nearest',fill_value="extrapolate")  # require to fit the reference phase to get H for any P values

for i in sys:
    if i == "FCC":
        print("\nReference phase is processing: "+i)
    else:
        print("Other phase is processing: "+i)
        extract = data(i,refint)
        datasave = np.append(datasave,extract,axis=0)
        
    
np.savetxt(f"{script_dir}/{xcfunc}/GlobalRel_H-{xcfunc}.dat",np.column_stack(datasave))  # all data in the same file according to this order : XI(P,H), II(P,H), XV(P,H), VIII(P,H)

####################################################################
# Searching for transition pressure

transitions = []

for i in sys:
    if i == "FCC":
        j = 0  # column for P values of ice phase XI
    else:
        j = j + 2
        refprev = interpolate.interp1d(
            datasave[j-2], datasave[j-1],
            kind='nearest', fill_value="extrapolate")
        P = datasave[j]
        H = datasave[j+1]
        order = np.argsort(P)
        P = P[order]
        H = H[order]
        diff = H - refprev(P)
        index = np.where(np.diff(np.sign(diff)))[0][0]
        ptr = P[index]
        transitions.append((ptr, sys[int(j/2-1)], i))
        print(
            f'Transition pressure between {sys[int(j/2-1)]} and {i}: '
            f'{ptr:.3f} GPa'
        )
      

###############################################################################
#Plottint the graph with all relative enthalpies
fig, ax = plt.subplots()
ax.plot(datasave[0],datasave[1],color='r',label='FCC')
ax.plot(datasave[2][::2],datasave[3][::2],color='b',label='HCP')
ax.set_xlabel('Pression (GPa)', fontweight='bold')
ax.set_ylabel('Relative enthalpy (eV)/molecule',fontweight='bold')
ax.set_xscale('log')
ax.set_xlim(0.01, 10)  
ax.set_ylim(-.005, .015) 

ax.minorticks_on()
leg =ax.legend(
    fontsize=12, 
    loc='upper left',  
    handlelength=1,  
    handletextpad=0.5,  
    borderpad=0.2,  
    labelspacing=0.2,  
    markerscale=0.8  
)
for text in leg.get_texts():
    text.set_fontweight('bold')
plt.tight_layout()

# Plot transition pressures on the graph

for ptr, phase1, phase2 in transitions:
    ax.axvline(
        x=ptr,
        color='k',
        linestyle='--',
        linewidth=1
    )
    
    ax.text(
        ptr*0.03, 0.0070,         
        f'{phase1} → {phase2}\n$P_{{tr}}$ = {ptr:.2f} GPa',
        rotation=0,
        fontsize=12,
        ha='left',
        va='top'
    )

plt.savefig(f"{script_dir}/{xcfunc}/relativeH-{xcfunc}.pdf")
plt.show()
