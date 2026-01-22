# data_DFT-thermodynamical-study-of-xenon-clathrate

Data from the article: "DFT Investigation of Xenon Clathrate Hydrates: thermodynamic properties and structural insights"

Here are some informations on what each folder contains.


* Convexhull:
Contains energy–volume data of xenon clathrate hydrate for various occupation rates, as well as for different pure reference phases (ICE + Xe). These data are needed for convex hull analysis. The attached numbers represent the following cases: 8: full occupation; 7-: under-occupation with one atom missing in the small cage; 7_: under-occupation with one atom missing in the large cage; 6-: under-occupation with two atoms missing in the small cage; and 6_: under-occupation with two atoms missing in the large cage. For pure reference phases, the data-XC-Phases.dat files contain information on the unit-cell parameter, volume, pressure, energy, enthalpy of the considered phases for a specific exchange-correlation (XC) functional. Python scripts to extract data-XC-Phases.dat files are also available, as well as scripts that compare the enthalpy between phases (relative enthalpy).


* Cell_parameter-temperature:
Evolution of the unit cell parameter with pressure for fully occupied xenon clathrate


* Thermal expansion and bulk modulus:
Contains data and scripts for analyzing the thermal expansion and bulk modulus of xenon clathrate hydrate. It includes: thermal expansion data from two reference articles, QHA data used for our calculations, and the Python scripts. Using these files, one can generate plots of the thermal expansion coefficient as a function of temperature and the bulk modulus B as a function of temperature for different pressures.
