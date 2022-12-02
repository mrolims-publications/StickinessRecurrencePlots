"""
Plot Figure 5. Simply run

    $ python plot_fig5.py

Make sure to run

    $ python fig5a.py

with Ntot = int(1e10) and

    $ python fig5bc_6a.py

with Ntot = int(1e8) and Ntot = int(1e9) before running this code.

Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import plot_params
import matplotlib as mpl
from scipy.interpolate import interp1d
import os

color1 = 'r'
color2 = 'lime'
color3 = 'gold'

color1 = 'dodgerblue'
color2 = 'red'
color3 = 'lime'

marker = 'x'
ms = 0.02

plot_params(fontsize=32, tick_labelsize=30, axes_labelsize=33)
mpl.rcParams['axes.linewidth'] = 1.6 #set the value globally
fig, ax = plt.subplots(1, 3, facecolor='w', figsize=(22, 6))

####################
# --- Fig 5(a) --- #
####################
k = 1.5
n = 200
Ntot = int(1e10)
exponent = int(np.log10(Ntot))
base = int(Ntot/10**exponent)
xbox = 0.0062
ybox = 0.93
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
path = 'Data/'
datafile = path + 'fig5a_Ntot=%ie%i.dat' % (base, exponent)
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig5a.py?\nStopping execution...' % (datafile))
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
X = np.array(df[0])

bin_heights, bin_borders, _ = ax[0].hist(X, density=True, histtype='step', color='k', bins='auto')
bin_widths = np.diff(bin_borders)
bin_centers = bin_borders[:-1] + bin_widths / 2
x = bin_centers
y = bin_heights
print("Interpolation the data...")
f = interp1d(x, y, kind='cubic')
print('Plotting the data o gráfico 1...')
s01a = 1.63
sf1a = 1.745
ax[0].fill_between(np.linspace(s01a, sf1a), f(np.linspace(s01a, sf1a)), alpha=0.95, color=color1)

s01b = sf1a
sf1b = 1.883
ax[0].fill_between(np.linspace(s01b, sf1b), f(np.linspace(s01b, sf1b)), alpha=0.95, color=color2)

s01c = sf1b
sf1c = 2.01
ax[0].fill_between(np.linspace(s01c, sf1c), f(np.linspace(s01c, sf1c)), alpha=0.95, color=color3)


s02 = 2.75
sf2 = 3.05
ax[0].fill_between(np.linspace(s02, sf2), f(np.linspace(s02, sf2)), alpha=0.95, color='k')
_ = ax[0].set_xlabel('$\\mathrm{RTE}(%i)$' % n), ax[0].set_ylabel('$P(\\mathrm{RTE}(%i))$' % n), ax[0].set_xlim(0, x.max()), ax[0].set_yticks([0, 0.5, 1, 1.5])
ax[0].text(xbox, ybox, '(a)', transform=ax[0].transAxes, bbox=bbox)
ax[0].set_xticks([0, 1, 2, 3, 4, 5])
# --- Inset --- #
ax_ins = ax[0].inset_axes([0.145, 0.4, 0.6, 0.45])
ax_ins.plot(np.arange(len(X)), X, 'k-', lw=0.1)
#ax[0].set_xscale('log')
ax_ins.set_xlim(40000, 70000)
ax_ins.set_xticks([40000, 70000])
#ax_ins.set_xticklabels(['$4.0\\times10^4$', '$8.0\\times10^4$'])
ax_ins.set_xlabel('$i$', fontsize=23)
ax_ins.set_ylabel('$\\mathrm{RTE}(200)$', fontsize=23)
ax_ins.set_yticks([0, 2.5, 5])
ax_ins.tick_params(axis='x', which='major', pad=10, labelsize=20)
ax_ins.tick_params(axis='y', which='major', labelsize=20)
print()

####################
# --- Fig 5(b) --- #
####################

k = 1.5
n = 200
Ntot = int(1e8)
exponent = np.log10(Ntot)
base = int(Ntot/10**exponent)

datafile1a = path + 'fig5b_blue_Ntot=%ie%i.dat' % (base, exponent)
datafile1b = path + 'fig5b_red_Ntot=%ie%i.dat' % (base, exponent)
datafile1c = path + 'fig5b_green_Ntot=%ie%i.dat' % (base, exponent)
datafile2 = path + 'fig5b_black_Ntot=%ie%i.dat' % (base, exponent)

# Checks if datafile exists
if not os.path.isfile(datafile2):
    import sys
    print('%s does not exist!\nDid you run fig5bc_6a.py?\nStopping execution...' % (datafile2))
    sys.exit()
print('Extracting data from %s...' % datafile2)
df = pd.read_csv(datafile2, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
print('Plotting the data...')
ax[1].plot(x, y, 'x', c='k', markersize=ms)

print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile1c, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
print('Plotting the data...')
ax[1].plot(x, y, 'x', c=color3, markersize=ms)

print('Extracting data from %s...' % datafile1b)
df = pd.read_csv(datafile1b, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
print('Plotting the data...')
ax[1].plot(x, y, 'x', c=color2, markersize=ms)

print('Extracting data from %s...' % datafile1a)
df = pd.read_csv(datafile1a, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
print('Plotting the data...')
ax[1].plot(x, y, 'x', c=color1, markersize=ms)


ax[1].text(xbox, ybox, '(b)', transform=ax[1].transAxes, bbox=bbox)

LW = 1.5
color = 'cyan'
ax[1].plot([1.45, 2.2], [-0.2, -0.2], '-', c=color, lw=LW)
ax[1].plot([2.2, 2.2], [-0.2, 1.7], '-', c=color, lw=LW)
ax[1].plot([1.45, 2.2], [1.7, 1.7], '-', c=color, lw=LW)
ax[1].plot([1.45, 1.45], [1.7, -0.2], '-', c=color, lw=LW)

_ = ax[1].set_xlim(-np.pi, np.pi), ax[1].set_ylim(-np.pi, np.pi)
_ = ax[1].set_xticks([-np.pi, 0, np.pi]), ax[1].set_yticks([-np.pi, 0, np.pi])
_ = ax[1].set_xticklabels(['$-\\pi$', '$0$', '$\\pi$']), ax[1].set_yticklabels(['$-\\pi$', '$0$', '$\\pi$']), 
_ = ax[1].set_ylabel('$p$')
_ = ax[1].set_xlabel('$x$')
print()

####################
# --- Fig 5(c) --- #
####################

ms = ms/2
k = 1.5
n = 200
Ntot = int(1e9)
exponent = np.log10(Ntot)
base = int(Ntot/10**exponent)

xi = 1.45
xf = 2.2
yi = -0.2
yf = 1.7

datafile1a = path + 'fig5b_blue_Ntot=%ie%i.dat' % (base, exponent)
datafile1b = path + 'fig5b_red_Ntot=%ie%i.dat' % (base, exponent)
datafile1c = path + 'fig5b_green_Ntot=%ie%i.dat' % (base, exponent)
datafile2 = path + 'fig5b_black_Ntot=%ie%i.dat' % (base, exponent)

# Checks if datafile exists
if not os.path.isfile(datafile2):
    import sys
    print('%s does not exist!\nDid you run fig5a.py?\nStopping execution...' % (datafile2))
    sys.exit()
print('Extracting data from %s...' % datafile2)
df = pd.read_csv(datafile2, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
X = x[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
Y = y[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
x = X
y = Y
print('Plotting the data...')
ax[2].plot(x, y, 'x', c='k', markersize=ms)

print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile1c, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
X = x[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
Y = y[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
x = X
y = Y
print('Plotting the data...')
ax[2].plot(x, y, 'x', c=color3, markersize=ms)

print('Extracting data from %s...' % datafile1b)
df = pd.read_csv(datafile1b, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
X = x[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
Y = y[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
x = X
y = Y
print('Plotting the data...')
ax[2].plot(x, y, 'x', c=color2, markersize=ms)

print('Extracting data from %s...' % datafile1a)
df = pd.read_csv(datafile1a, header=None, delim_whitespace=True)
x = np.array(df[0])
y = np.array(df[1])
X = x[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
Y = y[np.where((x >= xi) & (x <= xf) & (y >= yi) & (y <= yf))]
x = X
y = Y
print('Plotting the data...')
ax[2].plot(x, y, 'x', c=color1, markersize=ms)


ax[2].text(xbox, ybox, '(c)', transform=ax[2].transAxes, bbox=bbox)

_ = ax[2].set_xlim(xi, xf), ax[2].set_ylim(yi, yf)
_ = ax[2].set_xticks([xi, xf]), ax[2].set_yticks([yi, yf])
_ = ax[2].set_xlabel('$x$')

figname = 'Figures/fig5.png'
ax[0].tick_params(axis='x', which='major', pad=10)
ax[1].tick_params(axis='x', which='major', pad=10)
ax[2].tick_params(axis='x', which='major', pad=10)
_ = plt.subplots_adjust(left=0.05, bottom=0.18, right=0.98, top=0.966, hspace=0.27, wspace=0.2)
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=300, format='png')
print('Done.')
