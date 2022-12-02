"""
Generates the data of Figure 4(f). Simply run

    $ python fig4f.py

Execution time: hours
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from joblib import Parallel, delayed # Module to create parallel loops
from functions import RTE_border
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
# Creates the grid of points (k, y)
k = np.linspace(k0, k1, L, endpoint=True)
p = np.linspace(p0, p1, L, endpoint=True)
k, p = np.meshgrid(k, p)
# Computes the RTE for the grid above
rte = Parallel(n_jobs=-1)(delayed(RTE_border)(x, p[i, j], k[i, j], T) for i in range(L) for j in range(L))
rte = np.array(rte).reshape((L, L))
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig4f.dat'
# Saves the data
with open(datafile, 'w') as df:
    for i in range(L):
        for j in range(L):
            df.write('%.16f %.16f %.16f\n' % (k[i, j], p[i, j], rte[i, j]))
        df.write('\n')
