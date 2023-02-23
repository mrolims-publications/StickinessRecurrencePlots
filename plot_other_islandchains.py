import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import plot_params
import matplotlib as mpl
import os

color1 = 'r'
color2 = 'lime'
color3 = 'gold'

color1 = 'dodgerblue'
color2 = 'red'
color3 = 'lime'

marker = 'x'
ms = 0.01

ms = ms/2
n = 200
Ntot = int(1e9)
exponent = np.log10(Ntot)
base = int(Ntot/10**exponent)

path = '/home/matheus/Doutorado/StickinessRecurrencePlots/Data/'

datafile1a = path + 'fig5b_blue_Ntot=%ie%i.dat' % (base, exponent)
datafile1b = path + 'fig5b_red_Ntot=%ie%i.dat' % (base, exponent)
datafile1c = path + 'fig5b_green_Ntot=%ie%i.dat' % (base, exponent)
datafile2 = path + 'fig5b_black_Ntot=%ie%i.dat' % (base, exponent)


print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile1c, header=None, delim_whitespace=True)
x1 = np.array(df[0])
y1 = np.array(df[1])

print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile1b, header=None, delim_whitespace=True)
x2 = np.array(df[0])
y2 = np.array(df[1])

print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile1a, header=None, delim_whitespace=True)
x3 = np.array(df[0])
y3 = np.array(df[1])

print('Extracting data from %s...' % datafile1c)
df = pd.read_csv(datafile2, header=None, delim_whitespace=True)
x4 = np.array(df[0])
y4 = np.array(df[1])

xis = [1.45, 2.25, 2.1]
xfs = [2.2, np.pi, 2.7]
yis = [-0.2, -np.pi, -1.8]
yfs = [1.7, -1.9, -1.15]

plot_params(fontsize=32, tick_labelsize=30, axes_labelsize=33)
mpl.rcParams['axes.linewidth'] = 1.6 #set the value globally
fig, ax = plt.subplots(1, 3, facecolor='w', figsize=(22, 6))

label = ["(a)", "(b)", "(c)"]

xbox = 0.0062
ybox = 0.93
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
# Checks if datafile exists
if not os.path.isfile(datafile2):
    import sys
    print('%s does not exist!\nDid you run fig5a.py?\nStopping execution...' % (datafile2))
    sys.exit()

for i in range(len(xis)):

    xi = xis[i]
    xf = xfs[i]
    yi = yis[i]
    yf = yfs[i]

    
    X = x1[np.where((x1 >= xi) & (x1 <= xf) & (y1 >= yi) & (y1 <= yf))]
    Y = y1[np.where((x1 >= xi) & (x1 <= xf) & (y1 >= yi) & (y1 <= yf))]
    x = X
    y = Y
    ax[i].plot(x, y, 'x', c=color3, markersize=ms)

    X = x2[np.where((x2 >= xi) & (x2 <= xf) & (y2 >= yi) & (y2 <= yf))]
    Y = y2[np.where((x2 >= xi) & (x2 <= xf) & (y2 >= yi) & (y2 <= yf))]
    x = X
    y = Y
    ax[i].plot(x, y, 'x', c=color2, markersize=ms)

    X = x3[np.where((x3 >= xi) & (x3 <= xf) & (y3 >= yi) & (y3 <= yf))]
    Y = y3[np.where((x3 >= xi) & (x3 <= xf) & (y3 >= yi) & (y3 <= yf))]
    x = X
    y = Y
    ax[i].plot(x, y, 'x', c=color1, markersize=ms)

    X = x4[np.where((x4 >= xi) & (x4 <= xf) & (y4 >= yi) & (y4 <= yf))]
    Y = y4[np.where((x4 >= xi) & (x4 <= xf) & (y4 >= yi) & (y4 <= yf))]
    x = X
    y = Y
    ax[i].plot(x, y, 'x', c='k', markersize=ms/2)

    ax[i].text(xbox, ybox, label[i], transform=ax[i].transAxes, bbox=bbox)

    _ = ax[i].set_xlim(xi, xf), ax[i].set_ylim(yi, yf)
    _ = ax[i].set_xticks([xi, xf]), ax[i].set_yticks([yi, yf])
    _ = ax[i].set_xlabel('$x$')
    _ = ax[i].tick_params(axis='x', which='major', pad=10)


ax[0].set_ylabel("$p$")
ax[1].set_yticklabels(["$-\\pi$", "$-1.9$"])
ax[1].set_xticklabels(["$2.25$", "$\\pi$"])
figname = 'Figures/island_chains.png'
_ = plt.subplots_adjust(left=0.06, bottom=0.18, right=0.98, top=0.966, hspace=0.27, wspace=0.2)
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=250, format='png')