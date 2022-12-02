"""
Generates the data of Figure 4(a) and 4(b). Simply run

    $ python fig4ab.py a

to obtain the data from Figure 4(a) or

    $ python fig4ab.py b

to obtain the data from Figure 4(b).

Execution time: minutes
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from functions import lyapunov # Module with the standard map functions
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
    print('You must inform a or b.\ne.g. $ python fig4ab.py a')
    sys.exit()
idntf = sys.argv[1]
# Limits of the phase space
if idntf == 'a':
    x0 = -np.pi
    x1 = np.pi
    p0 = -np.pi
    p1 = np.pi
elif idntf == 'b':
    x0 = 1.45
    x1 = 2.2
    p0 = -0.2
    p1 = 1.7
else:
    print('Invalid identification. You must inform a or b.')
    sys.exit()
# Creates the grid of initial conditions
x = np.linspace(x0, x1, L, endpoint=True)
p = np.linspace(p0, p1, L, endpoint=True)
x, p = np.meshgrid(x, p)
# Evaluates the largest lyapunov exponent for the grid (x, p)
lyap = lyapunov(x, p, k, T)
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig4%s.dat' % (idntf)
# Saves the data
with open(datafile, 'w') as df:
    for i in range(L):
        for j in range(L):
            df.write('%.16f %.16f %.16f\n' % (x[i, j], p[i, j], lyap[i, j]))
        df.write('\n') # Add a blank line to plot in gnuplot using pm3d map