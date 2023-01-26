"""
Plot Figure 2. Simply run

    $ python plot_fig2.py

Make sure to run

    $ python fig2.py

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

# Define the parameters of the plot
plot_params(fontsize=21, tick_labelsize=23, axes_labelsize=26)
# (x, y) position of the label box
xbox = 0.0065
ybox = 0.9445
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
# Labels
label = ['(a)', '(b)', '(c)']
# Colors of the plot
color = ['b', 'k', 'r']
# Create the axes
fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, facecolor='w', figsize=(15, 5))
# Datafile
datafile = 'Data/fig2.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig2.py?\nStopping execution...' % datafile)
    sys.exit()
# Extract the data
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
# Plot the data
print('Plotting the data...')
for i in range(len(color)):
    recmat = np.array(df[i + 2])
    N = int(np.sqrt(len(recmat)))
    recmat = recmat.reshape((N, N))
    x = np.where(recmat == 1)[0]
    y = np.where(recmat == 1)[1]
    ax[i].scatter(x, y, c=color[i], marker='s', s=0.5, edgecolors='none')
    ax[i].text(xbox, ybox, label[i], transform=ax[i].transAxes, bbox=bbox)
    ax[i].set_xlabel('$i$')
# Adjust the plot
print('Adjusting the plot...')
mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally
ax[0].set_ylabel('$j$')
ax[0].set_xticks([0, 250, 500, 750, 1000])
ax[0].set_yticks([0, 250, 500, 750, 1000])
ax[0].set_xlim(0, N)
ax[0].set_ylim(0, N)
_ = plt.subplots_adjust(left=0.075, bottom=0.15, right=0.975, top=0.97, wspace=0.12)
# Save the figure in path
path = 'Figures/'
# Check if path exists
if not os.path.exists(path):
    os.system('mkdir %s' % path)
figname = path + 'fig2.png'
# Save the figure
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=300, format='png')
print('Done.')