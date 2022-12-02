"""
Plot Figure 6. Simply run

    $ python plot_fig6.py

Make sure to run

    $ python fig6b.py

with Ntot = int(1e11) and

    $ python fig5bc_6a.py

with Ntot = int(1e8) before running this code.

Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import plot_params
import matplotlib as mpl
import os

color1 = 'dodgerblue'
color2 = 'red'
color3 = 'lime'

plot_params(fontsize=27, tick_labelsize=25, axes_labelsize=28)
mpl.rcParams['axes.linewidth'] = 1.6 #set the value globally
fig, ax = plt.subplots(2, 1, facecolor='w', figsize=(8, 12))
xbox = 0.005
ybox = 0.94
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}

####################
# --- Fig 6(a) --- #
####################

Ntot = int(1e8)
exponent = np.log10(Ntot)
base = int(Ntot/10**exponent)
path = 'Data/'

datafile = path + 'fig6a_black_Ntot=%ie%i.dat' % (base, exponent)
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig5bc_6a.py?\nStopping execution...' % (datafile))
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
print('Plotting the data...')
ax[0].plot(x, y, 'kx', markersize=0.02)

ax[0].text(xbox, ybox, '(b)', transform=ax[1].transAxes, bbox=bbox)


_ = ax[0].set_xlim(-np.pi, np.pi), ax[0].set_ylim(-np.pi, np.pi)
_ = ax[0].set_xticks([-np.pi, 0, np.pi]), ax[0].set_yticks([-np.pi, 0, np.pi])
_ = ax[0].set_xticklabels(['$-\\pi$', '$0$', '$\\pi$']), ax[0].set_yticklabels(['$-\\pi$', '$0$', '$\\pi$']), 
_ = ax[0].set_ylabel('$p$')
_ = ax[0].set_xlabel('$x$')
ax[0].text(xbox, ybox, '(a)', transform=ax[0].transAxes, bbox=bbox)

####################
# --- Fig 6(b) --- #
####################

Ntot = int(1e11)
exponent = int(np.log10(Ntot))
base = int(Ntot/10**exponent)
datafile = path + 'fig6b_Ntot=%ie%i.dat' % (base, exponent)
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig6b.py?\nStopping execution...' % (datafile))
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
t1a = np.array(df[0])
Q1a = np.array(df[1])
t1b = np.array(df[2])
Q1b = np.array(df[3])
t1c = np.array(df[4])
Q1c = np.array(df[5])
t2 = np.array(df[6])
Q2 = np.array(df[7])
t3 = np.array(df[8])
Q3 = np.array(df[9])
ms = 2
print('Plotting the data...')
ax[1].plot(t1a, Q1a, 'o', c=color1, markersize=ms)
ax[1].plot(t1b, Q1b, 'o', c=color2, markersize=ms)
ax[1].plot(t1c, Q1c, 'o', c=color3, markersize=ms)
ax[1].plot(t3, Q3, 'k-', markersize=ms)
ax[1].plot(t2, Q2, 'ko', markersize=ms)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim(1e0, 1e4)
ax[1].set_xlabel('$\\tau$')
ax[1].set_ylabel('$Q(\\tau)$')
ax[1].text(xbox, ybox, '(b)', transform=ax[1].transAxes, bbox=bbox)
ax[1].set_ylim(1e-6, 1e0)
ax[1].set_xticks([1e0, 1e1, 1e2, 1e3, 1e4])
ax[1].set_yticks([1e0, 1e-2, 1e-4, 1e-6])

_ = plt.subplots_adjust(left=0.15, bottom=0.07, right=0.96, top=0.975, hspace=0.22, wspace=0.25)
figname = 'Figures/fig6.png'
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=300, format='png')
print('Done.')
