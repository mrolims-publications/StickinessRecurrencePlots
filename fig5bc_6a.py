"""
Generates the data of Figures 5(b), 5(c) and 6(a). Simply run

    $ python fig5bc_6a.py 1e8
    $ python fig5bc_6a.py 1e9

Execution time: minutes
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from functions import RTE_border, stdmap
import os # Module to check if the directory exists
import sys

if len(sys.argv == 1):
    print('You must inform Ntot.\nUse 1e8 AND 1e9.')
    sys.exit()

# Non-linearity parameter
k = 1.5
# Finite-time
n = 200
# Total number of iterations
Ntot = int(float(sys.argv[1]))
exponent = int(np.log10(Ntot))
base = int(Ntot/10**exponent)
# Number of windows of size n
N = round(Ntot/n)
# Initial condition (it can be any in the chaotic sea)
x0 = -3.0
p0 = 0
# Fig 5bc - blue
s01a = 1.63
sf1a = 1.745
# Fig 5bc - red
s01b = sf1a
sf1b = 1.883
# Fig 5bc - green
s01c = sf1b
sf1c = 2.01
# Fig 5bc - black
s02 = 2.75
sf2 = 3.05
# Fig 6a - black
s03 = 4
sf3 = 4.9

# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the files
datafile1a = path + 'fig5b_blue_Ntot=%ie%i.dat' % (base, exponent)
datafile1b = path + 'fig5b_red_Ntot=%ie%i.dat' % (base, exponent)
datafile1c = path + 'fig5b_green_Ntot=%ie%i.dat' % (base, exponent)
datafile2 = path + 'fig5b_black_Ntot=%ie%i.dat' % (base, exponent)
datafile3 = path + 'fig6a_black_Ntot=%ie%i.dat' % (base, exponent)
# Opens the files
df1a = open(datafile1a, 'w')
df1b = open(datafile1b, 'w')
df1c = open(datafile1c, 'w')
df2 = open(datafile2, 'w')
df3 = open(datafile3, 'w')
# Sets the initial condtion and iterates the orbit
x = x0
p = p0
for i in range(N):
    # Gets the phase-space positions
    ts = stdmap(x, p, k, n)
    X = ts[:, 0]
    P = ts[:, 1]
    # Evaluate the RTE for the (x, p) IC
    ftrte, x, p = RTE_border(x, p, k, n, return_last_pos=True)
    # Checks if RTE is on the intervals defined above. If so,
    # writes its phase-space positions to the corresponding file
    if ftrte >= s01a and ftrte <= sf1a:
        for j in range(n):
            df1a.write('%.16f %.16f\n' % (X[j], P[j]))
    elif ftrte > s01b and ftrte <= sf1b:
        for j in range(n):
            df1b.write('%.16f %.16f\n' % (X[j], P[j]))
    elif ftrte > s01c and ftrte <= sf1c:
        for j in range(n):
            df1c.write('%.16f %.16f\n' % (X[j], P[j]))
    elif ftrte >= s02 and ftrte <= sf2:
        for j in range(n):
            df2.write('%.16f %.16f\n' % (X[j], P[j]))
    elif ftrte >= s03 and ftrte <= sf3 and Ntot == 1e8:
        # Does not save when Ntot = 1e9 due to high number of points
        # Uses Ntot = 1e8 to plot the phase space
        for j in range(n):
            df3.write('%.16f %.16f\n' % (X[j], P[j]))
# Closes the files
df1a.close()
df1b.close()
df1c.close()
df2.close()
df3.close()