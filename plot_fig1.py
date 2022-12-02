"""
Plot Figure 1. Simply run

    $ python plot_fig1.py

Make sure to run

    $ python fig1.py

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
plot_params(fontsize=21, tick_labelsize=19, axes_labelsize=22)
# (x, y) position of the label box
xbox = 0.0065
ybox = 0.9445
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
# Labels
label = ['(a)', '(b)', '(c)']
# Colors of the plot
color = ['b', 'k', 'r']
# Create the axes
fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, facecolor='w', figsize=(6, 5))
# Datafile
datafile = 'Data/fig1.dat'
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig1.py?\nStopping execution...' % datafile)
    sys.exit()
# Extract the data
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
# Plot the data
print('Plotting the data...')
for i in range(int(len(df.columns)/2)):
    x = np.array(df[2*i])
    y = np.array(df[2*i + 1])
    ax.plot(x, y, 'x', color=color[i], markersize=0.7)
# Adjust the plot
print('Adjusting the plot...')
ax.set_xlabel('$x$')
_ = ax.set_xlim(-np.pi, np.pi), ax.set_ylim(-np.pi, np.pi)
_ = ax.set_xticks([-np.pi, 0, np.pi]), ax.set_yticks([-np.pi, 0, np.pi])
_ = ax.set_xticklabels(['$-\\pi$', '$0$', '$\\pi$']), ax.set_yticklabels(['$-\\pi$', '$0$', '$\\pi$']), 
_ = ax.set_ylabel('$p$')
_ = plt.subplots_adjust(left=0.12, bottom=0.15, right=0.98, top=0.98, wspace=0.075)
mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally
# Save the figure in path
path = 'Figures/'
# Check if path exists
if not os.path.exists(path):
    os.system('mkdir %s' % path)
figname = 'Figures/fig1.png'
# Save the figure
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=300, format='png')
print('Done.')