"""
Generates the data of Figure 1. Simply run

    $ python fig1.py

Execution time: seconds
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""
 
from functions import stdmap, lyapunov # Module with the standard map functions
import os # Module to check if the directory exists

# Nonlinearity parameter
k = 1.5
# Length of the orbits
N = 85000
# Initial conditons (IC)
x0 = [1, 2.9, 1.6]
p0 = [0, 0, 0]
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig1.dat'
# Empty list to save the xy time series of each IC
xy = []
# Iterates the initial contidions and evaluate the
# largest Lyapunov exponent of each orbit
for i in range(len(x0)):
    xy.append(stdmap(x0[i], p0[i], k, N))
    lyap = lyapunov(x0[i], p0[i], k, N)
    print('lyap_%i = %.5f' % (i+1, lyap))
# Saves the data
with open(datafile, 'w') as df:
    for i in range(len(xy[0])):
        df.write('%.16f %.16f %.16f %.16f %.16f %.16f\n' % (xy[0][i, 0], xy[0][i, 1], xy[1][i, 0], xy[1][i, 1], xy[2][i, 0], xy[2][i, 1]))