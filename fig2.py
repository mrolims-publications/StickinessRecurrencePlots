"""
Generates the data of Figure 2. Simply run

    $ python fig2.py

Execution time: seconds
Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

from functions import stdmap # Module with the standard map functions
from pyunicorn.timeseries import RecurrencePlot as RP # Pyunicorn module to create the recurrence plots
import os # Module to check if the directory exists

# Non-linearity parameter
k = 1.5
# Length of the orbits
N = 1000
# Initial conditions (ICs)
x0 = [1, 2.9, 1.6]
p0 = [0, 0, 0]
# Path to where the data will be stored
path = 'Data/'
# Checks if path exists
if not os.path.exists(path):
    # If not, creates path
    os.system('mkdir %s' % path)
# Name of the file
datafile = path + 'fig2.dat'
# Empty list to store the recurrence matrices of each orbit
recmats = []
# Iterates the initial conditions and store the orbit's recurrence matrices
for i in range(len(x0)):
    time_series = stdmap(x0[i], p0[i], k, N)
    rp = RP(time_series, metric='supremum', normalize=False, threshold_std=10/100, silence_level=2)
    recmats.append(rp.recurrence_matrix())
    print(recmats[i].shape)
# Saves the data
with open(datafile, 'w') as df:
    for i in range(N + 1):
        for j in range(N + 1):
            df.write('%i %i %i %i %i\n' % (i, j, recmats[0][i, j], recmats[1][i, j], recmats[2][i, j]))
        df.write('\n')