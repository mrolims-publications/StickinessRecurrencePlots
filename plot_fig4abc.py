"""
Plot Figure 4(a)(b)(c). Simply run

    $ python plot_fig4abc.py

Make sure to run

    $ python fig4ab.py a
    $ python fig4ab.py b
    $ python fig4c.py

before running this code.

Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from functions import plot_params

cmap = plot_params(fontsize=21, tick_labelsize=23, axes_labelsize=26)
xbox = 0.0065
ybox = 0.9445
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
label = 'abc'
xl = ['x', 'x', 'k']
fig, ax = plt.subplots(1, 3, facecolor='w', figsize=(16, 5))
path = 'Data/'
for i in range(len(label)):
    datafile = path + 'fig4%s.dat' % label[i]
    # Checks if datafile exists
    if not os.path.isfile(datafile):
        import sys
        if label[i] == 'c':
            script = 'fig4c.py'
        else:
            script = 'fig4ab.py'
        print('%s does not exist!\nDid you run %s?\nStopping execution...' % (datafile, script))
        sys.exit()
    print('Extracting data from %s...' % datafile)
    df = pd.read_csv(datafile, header=None, delim_whitespace=True)
    x = np.array(df[0])
    y = np.array(df[1])
    z = np.array(df[2])
    M = int(np.sqrt(len(x)))
    x = x.reshape((M, M))
    y = y.reshape((M, M))
    z = z.reshape((M, M))
    print('Plotting the data...')
    hm = ax[i].pcolor(x, y, z, vmin=0, vmax=0.5, shading='auto', cmap=cmap)
    ax[i].set_xlim(x.min(), x.max())
    ax[i].set_ylim(y.min(), y.max())
    ax[i].set_xlabel('$%s$' % xl[i])
    ax[i].text(xbox, ybox, '(%s)' % label[i], transform=ax[i].transAxes, bbox=bbox)
    if i == 0:
        ax[i].set_xticks([-np.pi, 0, np.pi]), ax[i].set_yticks([-np.pi, 0, np.pi])
        ax[i].set_xticklabels(['$-\\pi$', '$0$', '$\\pi$']), ax[i].set_yticklabels(['$-\\pi$', '$0$', '$\\pi$'])
    elif i == 1:
        ax[i].set_xticks([x.min(), x.max()]), ax[i].set_yticks([y.min(), y.max()])
    else:
        ax[i].set_xticks([0, 1, 2, 3, 4, 5]), ax[i].set_yticks([-np.pi, 0, np.pi])
        ax[i].set_yticklabels(['$-\\pi$', '$0$', '$\\pi$'])
    ax[i].tick_params(axis='x', which='major', pad=10)
    if i == 0:
        LW = 1
        ax[i].plot([1.45, 2.2], [-0.2, -0.2], 'w-', lw=LW)
        ax[i].plot([2.2, 2.2], [-0.2, 1.7], 'w-', lw=LW)
        ax[i].plot([1.45, 2.2], [1.7, 1.7], 'w-', lw=LW)
        ax[i].plot([1.45, 1.45], [1.7, -0.2], 'w-', lw=LW)
ax[0].set_ylabel('$p$')

cbar_ax = fig.add_axes([0.935, 0.165, 0.008, 0.97-0.165])
cbar = fig.colorbar(hm, cax=cbar_ax, label='$\\lambda_{\\mathrm{max}}$')
plt.subplots_adjust(left=0.055, bottom=0.165, right=0.93, top=0.97, wspace=0.17)
ax[-1].plot([0, 5], [1.3, 1.3], 'w--', lw=0.75)
# Save the figure in path
path = 'Figures/'
# Check if path exists
if not os.path.exists(path):
    os.system('mkdir %s' % path)
figname = 'Figures/fig4abc.png'
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=250, format='png')
print('Done.')