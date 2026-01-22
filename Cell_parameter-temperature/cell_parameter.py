import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

# Style configuration for xmgrace-like plots
plt.style.use('default')
plt.rcParams.update({
    'figure.figsize': (6, 6),
    'font.family': 'Garamond',
    'font.size': 14,
    'axes.linewidth': 2,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 16,
    'axes.grid': True,
    'grid.linestyle': '--',
    'grid.alpha': 0.3,
    'lines.linewidth': 0.8,
    'lines.markersize': 6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    
})

script_dir = os.path.dirname(__file__)



#Analytic function form Ikeda
def a_th(T):
    return 11.833 + 4.9692e-5 * T + 1.7966e-6 * T**2

# Temperature values for the analytical curve
temp_courbe = np.arange(0, 301, 1)
a_courbe = a_th(temp_courbe)

#Upload data from files optained from QHA calculations and from litterature
temp_sim,a_sim= np.loadtxt(f"{script_dir}/exp_simu.txt", skiprows=1, unpack=True)
temp_exp2,a_exp2 = np.loadtxt(f"{script_dir}/ikeda_values.txt", skiprows=1, unpack=True)
temp_vdw,a_vdw=np.loadtxt(f"{script_dir}/vdWDF2/volume-temperature_81.dat", unpack=True)
temp_rev, a_rev = np.loadtxt(f"{script_dir}/revPBE/volume-temperature_81.dat", unpack=True)


# Plot data
fig, ax = plt.subplots()
line1=ax.plot(temp_courbe, a_courbe, color='r', linestyle='-',linewidth=1.5, label='Fitting polynomial [1]')
line2=ax.plot(temp_exp2, a_exp2, color='r', marker='H', linestyle='--', label='Ikeda et al[1]')
line3=ax.plot(temp_sim, a_sim, color='k', marker='D', linestyle='-.', label='Lattice dynamics [2]')
line4=ax.plot(temp_vdw, a_vdw, color='g', marker='s', linestyle='-', label='vdWDF2')
line5=ax.plot(temp_rev, a_rev, color='purple', marker='v', linestyle='-', label='revPBE-D3(0)')

ax.set_xlabel('Temperature (K)', fontweight='bold')
ax.set_ylabel('Unit cell parameter (Å)', fontweight='bold')
ax.set_xlim([0, 310])
ax.set_ylim([11.70, 12.25])
leg=ax.legend(
    handles=[line2[0], line3[0], line4[0], line5[0]],
    labels=[
        'Experimental',
        'Force field',
        'vdW-DF2',
        'revPBE-D3(0)'
    ],
    fontsize=14,
    loc='best',
    handlelength=1,
    handletextpad=0.5,
    borderpad=0.2,
    labelspacing=0.2,
    markerscale=0.5
)
for text in leg.get_texts():
    text.set_fontweight('bold')
plt.tight_layout()

plt.grid(which='major', linestyle='-', linewidth=0.7, color='grey', alpha=0.5)  
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig(f'{script_dir}/cell_parameter-temperature.pdf', dpi=300)
plt.show()

