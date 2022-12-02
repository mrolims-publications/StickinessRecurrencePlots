"""
Generates the data of Figure 4(d) and 4(e). Simply run

    $ python fig4ab.py d

to obtain the data from Figure 4(d) or

    $ python fig4ab.py e

to obtain the data from Figure 4(e).

Execution time: hours
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""


import numpy as np # NumPy module
from joblib import Parallel, delayed # Module to create parallel loops
from functions import RTE_border # Module with the standard map functions
import os # Module to check if the directory exists
import sys

# Nonlinearity parameter
k = 1.5
# Size of initial condition's grid L x L 
L = 2**10
# Length of the orbits
T = int(5e3)
# Checks if the user informed the identification for the data
if len(sys.argv == 1):
    print('You must inform d or e.\ne.g. $ python fig4de.py e')
    sys.exit()
idntf = sys.argv[1]
# Limits of the phase space
if idntf == 'd':
    x0 = -np.pi
    x1 = np.pi
    p0 = -np.pi
    p1 = np.pi
elif idntf == 'e':
    x0 = 1.45
    x1 = 2.2
    p0 = -0.2
    p1 = 1.7
else:
    print('Invalid identification. You must inform d or e.')
    sys.exit()
# Creates the grid of initial conditions
x = np.linspace(x0, x1, L, endpoint=True)
p = np.linspace(p0, p1, L, endpoint=True)
x, p = np.meshgrid(x, p)
# Evaluate the RTE
rte = Parallel(n_jobs=-1)(delayed(RTE_border)(x[i, j], p[i, j], k, T) for i in range(L) for j in range(L))
rte = np.array(rte).reshape((L, L))
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig4%s.dat' % idntf
# Saves the data
with open(datafile, 'w') as df:
    for i in range(L):
        for j in range(L):
            df.write('%.16f %.16f %.16f\n' % (x[i, j], p[i, j], rte[i, j]))
        df.write('\n')