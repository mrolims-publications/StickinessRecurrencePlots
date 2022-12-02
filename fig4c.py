"""
Generates the data of Figure 4(c). Simply run

    $ python fig4c.py

Execution time: minutes
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from functions import lyapunov
import os # Module to check if the directory exists

# x initial condition
x = 0.0
# Size of the grid L x L
L = 2**10
# Length of the orbits
T = int(5e3)
# Limits of the parameter space
p0 = -np.pi
p1 = np.pi
k0 = 0.0
k1 = 5.0
# Creates the grid of points (k, p)
k = np.linspace(k0, k1, L, endpoint=True)
p = np.linspace(p0, p1, L, endpoint=True)
k, p = np.meshgrid(k, p)
# Computes the largest Lyapunov exponent for the grid (k, p)
lyap = lyapunov(x, p, k, T)
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig4c.dat'
# Saves the data
with open(datafile, 'w') as df:
    for i in range(L):
        for j in range(L):
            df.write('%.16f %.16f %.16f\n' % (k[i, j], p[i, j], lyap[i, j]))
        df.write('\n')