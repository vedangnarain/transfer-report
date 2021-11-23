#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jun 30 13:51:22 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Convert simulation units (nM) to partial pressure (mmHg)

Tested in Python 3.7.4.

Parameters and methods sourced from: 
    
[1] D. R. Grimes et al., ‘Estimating oxygen distribution from vasculature in 
three-dimensional tumour tissue’, J. R. Soc. Interface., vol. 13, no. 116, p. 
20160070, Mar. 2014, doi: 10.1098/rsif.2016.0070.

[2] T. D. Lewin, P. K. Maini, E. G. Moros, H. Enderling, and H. M. Byrne, ‘The 
Evolution of Tumour Composition During Fractionated Radiotherapy: Implications 
for Outcome’, Bull Math Biol, vol. 80, no. 5, pp. 1207–1235, May 2018, doi: 
10.1007/s11538-018-0391-9.

Used to generate figure in Transfer Report.

"""

# Initialise libraries
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Define function to convert nM to mmHg
def convert_nM_to_mmHg(nM, K=22779):
    
    # Convert nM (moles of O2 per cubic decimetre) into kg per cubic metre
    concentration = nM*0.000000001*15.999
    
    # Calculate partial pressure (mmHg) using formula from [2]
    mmHg = K*concentration
    
    # Return partial pressure
    return mmHg
 
# Define function to convert mmHg to nM
def convert_mmHg_to_nM(mmHg, K=22779):
    
    # Convert partial pressure (mmHg) to concentration (kg per cubic metre)
    concentration = mmHg/K
    
    # Convert kg per cubic metre into nM
    nM = concentration/(0.000000001*15.999)
    
    # Return partial pressure
    return nM    

# Enter parameters from paper [1]
K = 22779
anoxic_pp = 0.8
hypoxic_pp = 10
tumour_boundary_pp = 100

# Set LaTex-style font
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 11})
#plt.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

# Plot concentration (nM) vs partial pressure (mmHg)
nM_range = np.linspace(0, 40000, 100)
pp_range = nM_range*0.000000001*15.999*K
fg = plt.figure()
#plt.title('Relationship between concentration and partial pressure (Grimes et al., 2014)')
plt.plot(nM_range, pp_range, c='black')
#plt.axhline(100, alpha=0.5, c='g', ls='--', label='tumour boundary')
#plt.axhline(hypoxic_pp, alpha=0.5, c='orange', ls='--', label='hypoxia')
#plt.axhline(anoxic_pp, alpha=0.5, c='red', ls='--', label='anoxia')
plt.fill_between(nM_range, 0, 0.8, facecolor='blue', alpha=0.75)
plt.fill_between(nM_range, 0.8, 10, facecolor='deepskyblue', alpha=0.75)
plt.fill_between(nM_range, 10, 15, facecolor='red', alpha=0.75)
plt.xlim(0, 40000)
plt.ylim(0,15)
plt.xlabel('concentration of oxygen (nM)')
plt.ylabel('partial pressure of oxygen (mmHg)')
#plt.legend()

import matplotlib.patches as mpatches

anoxia_patch = mpatches.Patch(color='blue', label='anoxia')
hypoxia_patch = mpatches.Patch(color='deepskyblue', label='hypoxia')
normoxia_patch = mpatches.Patch(color='red', label='normoxia')

plt.legend(handles=[normoxia_patch, hypoxia_patch, anoxia_patch, ])

plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/relationship between concentration and partial pressure.svg').expanduser()
fg.savefig(file_path, bbox_inches='tight', dpi=500)
file_path = Path('~/Desktop/Final Figures/relationship between concentration and partial pressure.png').expanduser()
fg.savefig(file_path, bbox_inches='tight', dpi=500)
