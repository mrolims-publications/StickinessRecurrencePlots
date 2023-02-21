"""
Plot Figure Appendix. Simply run

    $ python plot_fig_appendix.py

Make sure to run

    $ python fig_appendix.py

before running this code.

Author: Matheus Rolim Sales
Last modified: 21/02/2022
"""

import numpy as np # NumPy module
from functions import plot_params, stdmap
import os # Module to check if the directory exists
import matplotlib.pyplot as plt
import pandas as pd

# Datafile
datafile = "Data/cc_vs_eps.dat"
# Checks if datafile exists
if not os.path.isfile(datafile):
    import sys
    print('%s does not exist!\nDid you run fig_appendix.py?\nStopping execution...' % datafile)
    sys.exit()
# Compute the standard deviation
print("Calculating the standard deviations...")
# Number of points in k
L = 5000
# Length of the orbit
T = int(5e3)
# Initial condition (IC)
x = 0.0
p = 1.3
# Limits of k
k_ini = 0
k_end = 5
# Creates an array from k_ini to k_end (included) with L points
k = np.linspace(k_ini, k_end, L, endpoint=True)
std1 = np.zeros(L)
std2 = np.zeros(L)
std3 = np.zeros(L)
for i in range(len(k)):
    time_series = stdmap(x, p, k[i], T)
    std1[i] = np.std(time_series)
    std2[i] = max(np.std(time_series[:, 0]), np.std(time_series[:, 1]))
    std3[i] = np.sqrt(np.std(time_series[:, 0])**2 + np.std(time_series[:, 1])**2)
print('Extracting data from %s...' % datafile)
df = pd.read_csv(datafile, header=None, delim_whitespace=True)
x = np.array(df[0])
y1 = np.array(df[1])
y2 = np.array(df[2])
y3 = np.array(df[3])
print('Plotting the data...')
# Linewidth
lw = 0.7
# (x, y) position of the label box
xbox = 0.003
ybox = 0.9035
bbox = {'linewidth': 0.0, 'facecolor': 'white', 'alpha': 1.0, 'pad': 1}
# Define the parameters of the plot
plot_params(fontsize=26, tick_labelsize=28, axes_labelsize=31, legend_fontsize=20)
# Create the axes
fig, ax = plt.subplots(2, 1, figsize=(10, 8), facecolor="w")
ax[0].plot(k, std1, "k", lw=lw)
ax[0].plot(k, std2, "r", lw=lw)
ax[0].plot(k, std3, "b", lw=lw)
ax[0].set_xlabel("$k$")
ax[0].set_ylabel("$\\sigma$")
ax[0].set_xlim(k.min(), k.max())
ax[0].text(xbox, ybox, '(a)', transform=ax[0].transAxes, bbox=bbox)
ax[0].set_yticks([1, 1.5, 2, 2.5])

ax[1].plot(x*100, y1, "kx-", lw=lw, label="Concatenation approach")
ax[1].plot(x*100, y2, "rx-", lw=lw, label="Maximum norm")
ax[1].plot(x*100, y3, "bx-", lw=lw, label="Euclidean norm")
ax[1].set_xlabel("$\\epsilon$ (in units of %$\\sigma$)"), ax[1].set_ylabel("$\\rho_{\\lambda_{\\mathrm{max}}, \\mathrm{RTE}}$")
ax[1].set_xlim((x*100).min(), (x*100).max())
ax[1].set_xticks([1, 5, 10, 15])
ax[1].text(xbox, ybox, '(b)', transform=ax[1].transAxes, bbox=bbox)
ax[1].legend(loc="lower right", frameon=False)
# Adjust the plot
print('Adjusting the plot...')
plt.subplots_adjust(left=0.12, bottom=0.115, right=0.97, top=0.98, hspace=0.3)
# Save the figure in path
path = 'Figures/'
# Check if path exists
if not os.path.exists(path):
    os.system('mkdir %s' % path)
figname = path + "fig_appendix.png"
# Save the figure
print('Saving in %s...' % figname)
plt.savefig(figname, dpi=500, format="png")
print('Done.')