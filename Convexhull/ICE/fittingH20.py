### Data fitting ###

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


script_dir = os.path.dirname(__file__)

###############################################################################
# vol for tetragonal(XI): 1.73207*1.6277
# vol for hexagonal(II): 0.48188*np.sin(120)
# vol for rhombohedron(II): np.sqrt(1-3*np.cos(np.pi*113.1/180)**2+2*np.cos(np.pi*113.1/180)**3)
# vol for triclinic (XV): 1.0018*0.92908
# vol for tetragonal (VIII): 1.45

xcfunc = "revPBE"  # choix fonctionelle

sys = "Ih"  # choice of the ice phase !!!

if sys == "XI":
    coeffvol=1.73207*1.6277  # define the volume
    nbwater=8                # number of water molecules per unit cell
elif sys == "II":
    coeffvol=np.sqrt(1-3*np.cos(np.pi*113.1/180)**2+2*np.cos(np.pi*113.1/180)**3)
    nbwater=12
elif sys == "XV":
    coeffvol=1.0018*0.92908
    nbwater=10
elif sys == "VIII":
    coeffvol=1.45
    nbwater=8
elif sys == "Ih":
    coeffvol=1.414
    nbwater=4
###############################################################################
# Functions definition (Birch-Murnaghan EOS)
    
def part1(a,a0):
    return (a0/a)**2-1

def part2(a,a0):
    return 6-4*(a0/a)**2

def part3(a,a0):
    return (a0/a)**7-(a0/a)**5

def birchmurn(a,e0,a0,bm0,bm0d):
    return e0+9/16*a0**3*coeffvol*bm0*(part1(a,a0)**3*bm0d+part1(a,a0)**2*part2(a,a0))

def pressure(a,e0,a0,bm0,bm0d):
     return (3*bm0/2)*(part3(a,a0))*(1+(3/4*(bm0d-4)*(part1(a,a0))))
 
def enthalpy(a,e0,a0,bm0,bm0d):
    return (birchmurn(a,e0,a0,bm0,bm0d)+pressure(a,e0,a0,bm0,bm0d)*a**3*coeffvol)/nbwater


###############################################################################
# Fitting procedure using non-linear least squares

#e0=-832.9675
#a0=1727.4385
#bm0=6.8262
#bm0d=0.181201

data = np.loadtxt(script_dir+'/'+xcfunc+'/'+sys+'/energydft.dat', usecols=(0, 1))
adft = data[:,0]
edft = data[:,1]

popt,pcov = curve_fit(birchmurn,adft,edft,p0=[min(edft),min(adft),10,1])
perr = np.sqrt(np.diag(pcov))
xint = np.linspace(min(adft),max(adft),1000)
xpre = np.linspace(min(adft),popt[1],1000)


###############################################################################
# Saving data as : a,V,P,E,H

datasave=np.array([xpre,xpre**3*coeffvol/nbwater,pressure(xpre,*popt)*160.2,birchmurn(xpre,*popt)/nbwater,enthalpy(xpre,*popt)])
dataent=[xpre,enthalpy(xpre,*popt)]
np.savetxt(script_dir+'/'+xcfunc+'/'+sys+'/data-'+xcfunc+'-'+sys+'.dat',np.column_stack(datasave))
#np.savetxt(f"{script_dir}/{xcfunc}/{sys}/data-{xcfunc}-{sys}.dat",np.column_stack(datasave))


###############################################################################
# Fitting procedure using non-linear least squares

parameters = {'axes.labelsize': 7,'xtick.labelsize': 6,'ytick.labelsize': 6,'legend.fontsize': 5}
plt.rcParams.update(parameters)
plt.gcf().subplots_adjust(wspace=0.5)
plt.subplot(221)
plt.plot(adft, edft,marker="o",ls='',label="DFT results",markersize=4)
plt.plot(xint, birchmurn(xint,*popt),label="fitting curve")
#plt.title('Relaxation volume of H20-'+xcfunc, loc='left')
plt.ylabel('Potential energy (eV)')
plt.legend(loc="upper right")
plt.subplot(223)
plt.xlabel('Lattice parameter (Ang)')
plt.ylabel('Pressure (GPa)')
plt.plot(xpre,pressure(xpre,*popt)*160.2)

plt.subplot(222)
plt.ylabel('Potential energy (eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,birchmurn(xpre,*popt)/nbwater)
plt.subplot(224)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Enthalpy (eV/molec)')
plt.plot(pressure(xpre,*popt)*160.2,enthalpy(xpre,*popt))
#plt.suptitle('Relaxation volume of H20-'+xcfunc)
plt.savefig(script_dir+'/'+xcfunc+'/'+sys+'/properties-'+xcfunc+'-'+sys+'.pdf')
plt.show()


print(80*"-")
print("Fitted parameters are:",popt)
print("Error parameters are:", perr) 
print("\nVolume per H2O molecule:",popt[1]**3*coeffvol/nbwater,"+/-",3*perr[1]/popt[1]*popt[1]**3*coeffvol/nbwater,"Ang**3")
print("corresponding to a lattice parameter:", popt[1],"+/-", perr[1], "Ang" )
print("\nThe bulk modulus is:",popt[2]*160.2,"+/-",perr[2]*160.2,"GPa")


