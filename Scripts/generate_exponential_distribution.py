#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jun  1 18:23:28 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Generate a probability distribution function for P(death) as a function vessel diameters

Tested in Python 3.7.4.

Used to generate figure in Transfer Report.

"""

# Initialise libraries
import matplotlib.pyplot as plt
import numpy as np

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})
'''
x = np.linspace(0, 35, 100)
y1 = np.exp(-x/1)
y5 = np.exp(-x/5)
y10 = np.exp(-x/10)
y15 = np.exp(-x/15)
y20 = np.exp(-x/20)
y25 = np.exp(-x/25)
y30 = np.exp(-x/30)

fig = plt.figure()
plt.plot(x, y1, label='${γ}$ = 1', ls = 'solid')
plt.plot(x, y5, label='${γ}$ = 5', ls = 'dotted')
plt.plot(x, y10, label='${γ}$ = 10', ls = 'dashed')
#plt.plot(x, y15, label='${γ}$ = 15', ls = 'dashdot')
plt.plot(x, y20, label='${γ}$ = 20', ls= 'dashdot')
plt.plot(x, y30, label='${γ}$ = 30', ls = (0, (3, 1, 1, 1, 1, 1)))
plt.xlabel('vessel radius (μm)')
plt.ylabel('probability of death')
plt.legend()
plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/stochastic_pruning_probabilities.svg').expanduser()
fig.savefig(file_path, dpi=500)
file_path = Path('~/Desktop/Final Figures/stochastic_pruning_probabilities.png').expanduser()
fig.savefig(file_path, dpi=500)
'''
# Make plot for different sample vessel diameters
diameters = [1.25, 2.5, 3.5, 12.5, 30]
gammas = np.linspace(1, 26, 25)
y1 = np.exp(-diameters[0]/gammas)
y2 = np.exp(-diameters[1]/gammas)
y3 = np.exp(-diameters[2]/gammas)
y4 = np.exp(-diameters[3]/gammas)
y5 = np.exp(-diameters[4]/gammas)
#y30 = np.exp(-diameters[5]/gammas)

x = gammas

fig = plt.figure(figsize=(10,6), tight_layout = {'pad': 2})
plt.plot(x, y1, label='1.25 μm', ls = 'solid')
plt.plot(x, y2, label='2.5 μm', ls = 'dotted')
plt.plot(x, y3, label='3.5 μm', ls = 'dashed')
plt.plot(x, y4, label='12.5 μm', ls= 'dashdot')
plt.plot(x, y5, label='30 μm', ls = (0, (3, 1, 1, 1, 1, 1)))
plt.xlabel('${γ}$ (μm)')
plt.ylabel('probability of death')
#plt.xlim(1,25)
plt.legend(loc='best', prop={'size': 15})
plt.show()

# Save image
file_path = Path('~/Desktop/Transfer Material/Final Figures/vessel_stochastic_pruning_probabilities.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Transfer Material/Final Figures/vessel_stochastic_pruning_probabilities.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')

