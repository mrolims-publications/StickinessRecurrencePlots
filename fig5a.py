"""
To generate the data of Figure 5(a). Simply run

    $ python fig5a.py 1e10

To generate the data used in Figure 6(b), run

    $ python fig5a.py 1e11

Execution time: several hours for 1e10 and about 1.5 day for 1e11 (!!!!)
Author: Matheus Rolim Sales
Last modified: 05/12/2022
"""

import numpy as np # NumPy module
from functions import FTRTE_border
import os # Module to check if the directory exists
import sys

if len(sys.argv) != 3:
    print('You must inform Ntot and n.\nUse Ntot = 1e10 and n = 200 to generate the data from Figure 5(a) and Ntot = 1e11 and n = 200 to generate the data used in Figure 6(b).')
    sys.exit()
# Non-linearity parameter
k = 1.5
# Total number of iterations
Ntot = int(float(sys.argv[1]))
exponent = int(np.log10(Ntot))
base = int(Ntot/10**exponent)
# Finite-time
n = int(sys.argv[2]) #200
# Number of windows of size n
N = round(Ntot/n)
# Initial condition (it can be any in the chaotic sea)
x0 = -3.0
p0 = 0
# Evaluates the FTRTE distribuition
ftrte = FTRTE_border(x0, p0, k, n, Ntot)
# Path to where the data will be stored
path = 'Data/'
#path = '/home/matheus/StickinessRecurrencePlots/Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig5a_Ntot=%ie%i_n=%i.dat' % (base, exponent, n)
# Save the data
np.savetxt(datafile, ftrte)