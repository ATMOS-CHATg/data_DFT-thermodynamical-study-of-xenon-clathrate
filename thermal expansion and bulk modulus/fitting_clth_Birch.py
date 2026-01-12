#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:12:31 2024

@author: meko
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from temp import inject_temperature_block_with_transform
from scipy.optimize import curve_fit
import os

script_dir = os.path.dirname(__file__)

# Style général proche xmgrace
mpl.rcParams.update({
    "axes.linewidth": 2,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 14,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "font.family": "Garamond",
    "figure.dpi": 120,
    "figure.figsize": (8, 6),
    "savefig.dpi": 300,
})

# Load data from file
xcfunc = "vdWDF2"

source_file = f"{script_dir}/helmholtz-volume_{xcfunc}.dat"
target_temp_file = f"{script_dir}/{xcfunc}/energydft.dat"
temperature_list = np.arange(0, 301.0, 50).tolist()

all_pressures = []
all_bulk = []
all_xpre = []
all_T = []

for T in temperature_list:
    used_T = inject_temperature_block_with_transform(T, source_file, target_temp_file)
    data = np.loadtxt(f'{script_dir}/{xcfunc}/energydft.dat')
    a = data[:, 0]
    E_values = data[:, 1]

    # Define the Birch-Murnaghan equation
    def birch_murnaghan_fit(a,E0,a0,B0,dB0):
        term1=9*a0**3*B0/16
        term2=((a0/a)**2-1)**3
        term3=((a0/a)**2-1)**2
        term4=(6-4*(a0/a)**2)
        Energy = E0 + term1*(term2*dB0 + term3*term4)
        return Energy

    # Define pressure
    def pressure(a,E0,a0,B0,dB0):
        term1=((a0/a)**7-(a0/a)**5)
        term2=(3*(dB0-4)/4)
        term3=((a0/a)**2-1)
        return 3*B0/2*term1*(1+term2*(term3))

    #Define bulk modulus
    def bulk_modulus(a, a0, B0, dB0):
        term1 = (1/2) * B0 * (7 * (a / a0) ** (-7/3) - 5 * (a / a0) ** (-5/3))
        term2 = (3/8) * B0 * dB0 * (9 * (a / a0) ** (-3) - 14 * (a / a0) ** (-7/3) + 5 * (a / a0) ** (-5/3))
        Bulk = term1 + term2
        return Bulk

    popt, pcov = curve_fit(birch_murnaghan_fit, a, E_values, p0=[min(E_values), np.mean(a), 1, 0.1])
    perr = np.sqrt(np.diag(pcov))
    xpre = np.linspace(min(a), popt[1], 5000)

    pressure_values = pressure(xpre, popt[0], popt[1], popt[2]*160.2, popt[3])
    bulk_values = bulk_modulus(xpre, popt[1], popt[2]*160.2, popt[3])

    # Enregistre chaque résultat dans un fichier séparé
    np.savetxt(f'{script_dir}/{xcfunc}/PVB_results_T{int(T)}.dat',
               np.column_stack((xpre, pressure_values, bulk_values)),
               header='#Lattice_Parameter(Ang) Pressure(GPa) Bulk_Modulus(GPa)',
               fmt='%.6f',
               comments='')

    # Stocke pour le graphe final
    all_pressures.append(pressure_values)
    all_bulk.append(bulk_values)
    all_xpre.append(xpre)
    all_T.append(T)

    

# Après la boucle, trace tout sur un même graphe
plt.figure(figsize=(6, 6))
for xpre, pressure_values, T in zip(all_xpre, all_bulk, all_T):
    plt.plot(pressure_values, xpre, label=f'T = {int(T)} K')

plt.xlabel('Pressure (GPa)', fontweight='bold')
plt.ylabel('Unit cell parameter (Å)', fontweight='bold')
plt.grid(which='major', linestyle='-', linewidth=0.7, color='grey', alpha=0.5)

leg = plt.legend(
    fontsize=15,
    loc='best',
    handlelength=1,
    handletextpad=0.5,
    borderpad=0.2,
    labelspacing=0.2,
    markerscale=0.8
)
for text in leg.get_texts():
    text.set_fontweight('bold')
plt.tight_layout()

ax = plt.gca()
ax.margins(x=0, y=0)      # Supprime les marges
#ax.set_xlim(left=0)       # Commence l'axe X à 0
#ax.set_ylim(bottom=11.30)  # Décommente et adapte si tu veux aussi pour Y

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

plt.savefig(f"{script_dir}/{xcfunc}/Pressure_vs_Bulk_AllTemps.pdf")
plt.show()






