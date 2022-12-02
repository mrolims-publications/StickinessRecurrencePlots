"""
Plot Figure 3. Simply run

    $ python plot_fig3.py

Make sure to run

    $ python fig3.py

before running this code.

Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import os
from functions import plot_params

plot_params(fontsize=26, xtick_labelsize=28, ytick_labelsize=28, axes_labelsize=31)
xbox = 0.0065
ybox = 0.9
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
label = ['(a)', '(b)']

fig, ax = plt.subplots(2, 1, sharex=True, facecolor='w', figsize=(10, 7))
# Datafile
datafile = 'Data/fig3.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig3.py?\nStopping execution...' % datafile)
    sys.exit()
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
k = np.array(df[0])
lyap = np.array(df[1])
rte = np.array(df[2])
print('Plotting the data...')
ax[0].plot(k, lyap, 'k-', lw=0.7)
ax[1].plot(k, rte, 'k-', lw=0.7)
print('Adjusting the plot...')
ax[1].set_xlim(k.min(), k.max())
ax[0].text(xbox, ybox, label[0], transform=ax[0].transAxes, bbox=bbox)
ax[0].set_ylabel('$\\lambda$')
ax[1].text(xbox, ybox, label[1], transform=ax[1].transAxes, bbox=bbox)
ax[1].set_ylabel('RTE')
ax[1].set_xlabel('$k$')
mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally
_ = plt.subplots_adjust(left=0.105, bottom=0.12, right=0.98, top=0.97, hspace=0.05)
# Save the figure in path
path = 'Figures/'
# Check if path exists
if not os.path.exists(path):
    os.system('mkdir %s' % path)
figname = path + 'fig3.png'
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=300, format='png')
print('Evaluating the correlation coefficient...')
std_y1 = np.std(lyap)
std_y2 = np.std(rte)
covy1y2 = np.cov(lyap, rte)[0][1]
cc = covy1y2/(std_y1*std_y2)
print('Correlation coefficient: %.10f' % cc)
print('Done.')