"""
Generates the data of Figure 5(a). Simply run

    $ python fig5a.py

with Ntot = int(1e10) and Ntot = int(1e11).

Execution time: several hours (!!!!)
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from functions import FTRTE
import os # Module to check if the directory exists

# Non-linearity parameter
k = 1.5
# Finite-time
n = 200
# Total number of iterations
Ntot = int(1e11)
exponent = int(np.log10(Ntot))
base = int(Ntot/10**exponent)
# Number of windows of size n
N = round(Ntot/n)
# Initial condition (it can be any in the chaotic sea)
x0 = -3.0
p0 = 0
# Evaluates the FTRTE distribuition
ftrte = FTRTE(x0, p0, k, n, Ntot)
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig5a_Ntot=%ie%i.dat' % (base, exponent)
# Save the data
np.savetxt(datafile, ftrte)