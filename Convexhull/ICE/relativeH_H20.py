### Data fitting ###

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os

script_dir = os.path.dirname(__file__)
print (script_dir)

xcfunc = "vdWDF2" # choice of the XC functional

def data(phase,ref):
    a = np.loadtxt(script_dir+'/'+xcfunc+'/'+phase+'/data-'+xcfunc+'-'+phase+'.dat')

    return a[:,2],a[:,4]-ref(a[:,2])


sys=["XI","II","XV","VIII"]  # choice of the ice phase with ice XI as reference !!!


###############################################################################
# Calculation of relative enthalpies
print(f"{script_dir}/{xcfunc}/XI/data-{xcfunc}-XI.dat")
dataref = np.loadtxt(f"{script_dir}/{xcfunc}/XI/data-{xcfunc}-XI.dat")
dataref = dataref[:,2],dataref[:,4]
datasave = dataref[0],dataref[1]-dataref[1]
refint = interpolate.interp1d(dataref[0],dataref[1],kind='nearest',fill_value="extrapolate")  # require to fit the reference phase to get H for any P values

for i in sys:
    if i == "XI":
        print("\nReference phase is processing: "+i)
    else:
        print("Other phase is processing: "+i)
        extract = data(i,refint)
        datasave = np.append(datasave,extract,axis=0)
    
np.savetxt(f"{script_dir}/{xcfunc}/GlobalRel_H-{xcfunc}.dat",np.column_stack(datasave))  # all data in the same file according to this order : XI(P,H), II(P,H), XV(P,H), VIII(P,H)


###############################################################################
#Plottint the graph with all relative enthalpies

plt.plot(datasave[0],datasave[1],color='k',label='XI')
plt.plot(datasave[2],datasave[3],color='r',label='II')
plt.plot(datasave[4],datasave[5],color='g',label='XV')
plt.plot(datasave[6],datasave[7],color='b',label='VIII')
#plt.plot(datasave[8],datasave[9],color='y',label='Ih')
plt.xlabel('Pression(GPa) ')
plt.ylabel('ReEnthalpy(eV)')
plt.xscale('log')
plt.xlim(0, 6)  # Interval sur l'axe x de 0 à 6
plt.ylim(-.3, .08) 
plt.title(f"Relative enthalpy of ice phases with {xcfunc} functional")
plt.grid()

plt.legend()
plt.savefig(f"{script_dir}/{xcfunc}/relativeH-{xcfunc}_6kpts_600ev22.pdf")
plt.show()
###############################################################################
# Searching for transition pressure

print("\n"+"*"*50)
print("\n Searching for transition pressure:")
for i in sys:
    if i == "XI":
        j = 0  # column for P values of ice phase XI
    else:
        j = j+2
        refprev = interpolate.interp1d(datasave[j-2],datasave[j-1],kind='nearest',fill_value="extrapolate")
        diff = np.absolute(datasave[j+1]-refprev(datasave[j]))
        index = diff.argmin()
        ptr = datasave[j,index]
        print(f'Transition pressure between {sys[int(j/2-1)]} and {i}: {ptr:.3f} GPa')
#        
#        
#        
#
















