
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from temp import inject_temperature_block_with_transform
from scipy.optimize import curve_fit
from scipy import interpolate
from matplotlib.ticker import ScalarFormatter
import os

script_dir = os.path.dirname(__file__)


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
    "legend.fontsize": 12,
    "legend.loc": 'best',
    "legend.frameon": True,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "font.family": "Garamond",
    "figure.dpi": 120,
    "figure.figsize": (8, 6),
    "savefig.dpi": 300,
})


# Parameters
xcfunc = "revPBE"
source_file = f"{script_dir}/helmholtz-volume_{xcfunc}.dat"
target_temp_file = f"{script_dir}/{xcfunc}/energydft.dat"
temperature_list = np.arange(0, 301, 10).tolist() 

# Common pressure grid for interpolation
common_pressures = np.arange(0, 2.1, 0.25)  

# Dictionary for storing a(T) at each pressure
results_byP = {P: [] for P in common_pressures}
alpha_byP = {}

# Loop over temperature values
for T in temperature_list:
    used_T = inject_temperature_block_with_transform(T, source_file, target_temp_file)
    data = np.loadtxt(f'{script_dir}/{xcfunc}/energydft.dat')
    a = data[:, 0]
    E_values = data[:, 1]

    # Birch-Murnaghan
    def birch_murnaghan_fit(a, E0, a0, B0, dB0):
        term1 = 9*a0**3*B0/16
        term2 = ((a0/a)**2 - 1)**3
        term3 = ((a0/a)**2 - 1)**2
        term4 = (6 - 4*(a0/a)**2)
        return E0 + term1*(term2*dB0 + term3*term4)

    # Pression issue du fit
    def pressure(a, E0, a0, B0, dB0):
        term1 = ((a0/a)**7 - (a0/a)**5)
        term2 = (3*(dB0-4)/4)
        term3 = ((a0/a)**2 - 1)
        return 3*B0/2*term1*(1 + term2*term3)
    
    def bulk_modulus(a, a0, B0, dB0):
        term1 = (1/2) * B0 * (7 * (a / a0) ** (-7/3) - 5 * (a / a0) ** (-5/3))
        term2 = (3/8) * B0 * dB0 * (9 * (a / a0) ** (-3) - 14 * (a / a0) ** (-7/3) + 5 * (a / a0) ** (-5/3))
        Bulk = term1 + term2
        return Bulk

    # Fit
    popt, pcov = curve_fit(birch_murnaghan_fit, a, E_values, p0=[min(E_values), np.mean(a), 1, 0.1])
    xpre = np.linspace(min(a)-0.25, popt[1], 5000)
    perr=np.sqrt(np.diag(pcov))
    pressure_values = pressure(xpre, popt[0], popt[1], popt[2]*160.2, popt[3])
    bulk_values = bulk_modulus(xpre, popt[1], popt[2]*160.2, popt[3])

    # Interpolation on the common pressure grid
    f_interp = interpolate.interp1d(pressure_values, xpre, bounds_error=False, fill_value=np.nan)
    f_interp_B = interpolate.interp1d(pressure_values, bulk_values, bounds_error=False, fill_value=np.nan)
     
    for P in common_pressures:
        a_at_P = f_interp(P)
        B_at_P = f_interp_B(P)
        results_byP[P].append((T, a_at_P, B_at_P))

# Thermal expansion coefficient α(T) calculated for each pressure

def thermal_expansion(a_vals, T_vals):
    a0 = a_vals[0]
    alpha = np.zeros_like(a_vals)
    for i in range(1, len(T_vals)):
        alpha[i] = (a_vals[i] - a_vals[i-1]) / (a0 * (T_vals[i] - T_vals[i-1]))
    return alpha


for P, values in results_byP.items():
    T_vals = np.array([v[0] for v in values])
    a_vals = np.array([v[1] for v in values])
    B_vals = np.array([v[2] for v in values])
    a0 = a_vals[0]  # valeur à T=0
    alpha_vals = thermal_expansion(a_vals, T_vals)
    alpha_byP[P] = (T_vals,a_vals, B_vals,alpha_vals)

    
# alpha_ikeda
def a_th(T):
    return 11.833 + 4.9692e-5 * T + 1.7966e-6 * T**2

def alpha_th(T):
    da_dT = 4.9692e-5 + 2 * 1.7966e-6 * T
    return da_dT / a_th(T)

T,a = np.loadtxt(f"{script_dir}/ikeda_values.txt", skiprows=1, unpack=True)
V=a**3
da = np.diff(a)
dT = np.diff(T)
dadT = da / dT                 # a = cell parameter, T = temperature
a_i = a[:-1]
T_milieu = (T[:-1] + T[1:]) / 2

beta = dadT / a_i
print (beta)

 #  Import Hansen's experimental data
temp_exp,alpha_exp = np.loadtxt(f"{script_dir}/alpha_hansen.dat", skiprows=1, unpack=True)
alpha_exp2 =alpha_exp*1e-5


   # Plot thermal expansion coefficient α(T) for each pressure"
plt.figure(figsize=(6,6))
for P, (T_vals,a_vals, B_vals,alpha_vals) in alpha_byP.items():
    plt.plot(T_vals, alpha_vals, label=f"P={P:.2f} GPa")
plt.xlabel("Temperature (K)", fontweight='bold')
plt.xlim(0, 300)
plt.ylim(-0.6e-4,  1.20e-4)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.ylabel("Thermal expansion coefficient (1/K)", fontweight='bold')
plt.plot(temp_exp,alpha_exp2, color='k', marker='8', linestyle='', label='Hansen')
plt.plot(T_milieu,beta, color='m', marker='>', linestyle='', label='ikeda')

leg = plt.legend(
    ncol=3,
    frameon=True,              
    fancybox=True,                         
    framealpha=0.4,              
    edgecolor="black",
    handlelength=0.5,
    handletextpad=0.5,
    borderpad=0.5,
    labelspacing=0.05,
    markerscale=0.5
)
for text in leg.get_texts():
    text.set_fontweight('bold')
plt.tight_layout()

ax = plt.gca() 

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')
    
plt.grid(which='major', linestyle='-', linewidth=0.7, color='grey', alpha=0.5)
plt.savefig(f"{script_dir}/{xcfunc}/alpha_vs_Temp_Allpressure.pdf")
plt.show()

 #  Plot bulk modulus B(T) for each pressure
 
plt.figure(figsize=(6,6))
for P, (T_vals,a_vals, B_vals,alpha_vals) in alpha_byP.items():
    plt.plot(T_vals, B_vals, label=f"P={P:.2f} GPa")
    
plt.xlabel("Temperature (K)", fontweight='bold')
plt.ylabel("Bulk modulus (GPa)", fontweight='bold')
plt.xlim(0, 300)
plt.ylim(7, 15)

leg = plt.legend(
    ncol=2,
    handlelength=1,
    handletextpad=0.5,
    borderpad=0.2,
    labelspacing=0.05,
    markerscale=0.8
)
for text in leg.get_texts():
    text.set_fontweight('bold')
plt.tight_layout()

ax = plt.gca() 

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontweight('bold')
    
    
plt.grid(which='major', linestyle='-', linewidth=0.7, color='grey', alpha=0.5)
plt.savefig(f"{script_dir}/{xcfunc}/Bulk_vs_Temp_Allpressure.pdf")
plt.show()

# save data to a file
with open(f"{script_dir}/{xcfunc}/B-alpha_vs_T_for_eachP.dat", "w") as f:
    f.write("#Pressure(GPa) Temperature(K) Bulk_Modulus(GPa) Thermal_Expansion(1/K)\n")
    for P, (T_vals,a_vals, B_vals, alpha_vals) in alpha_byP.items():
        for i in range(len(T_vals)):
            f.write(f"{P:.6f} {T_vals[i]:.1f} {a_vals[i]:.6e} {B_vals[i]:.6f} {alpha_vals[i]:.6e}\n")

print(" Fichier généré :", f"{script_dir}/{xcfunc}/B-alpha_vs_T_for_eachP.dat")
