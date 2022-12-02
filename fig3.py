"""
Generates the data of Figure 3. Simply run

    $ python fig3.py


Execution time: minutes
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from joblib import Parallel, delayed # Module to create parallel loops
from functions import lyapunov, RTE # Module with the standard map functions
import os # Module to check if the directory exists

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
# Obtains the largest Lyapunov exponent for the array above using parallel computing with Numba.
lyap = lyapunov(x, p, k, T)
# Obtains the rte for the array above using parallel computing.
# n_jobs=-1 uses all available threads. Set it to another value if you 
# wish to use less.
rte = Parallel(n_jobs=-1)(delayed(RTE)(x, p, k[i], T) for i in range(L))
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig3.dat'
# Saves the data
with open(datafile, 'w') as df:
    for i in range(L):
        df.write('%.16f %.16f %.16f\n' % (k[i], lyap[i], rte[i]))