"""
Generates the data of Figure 7. Simply run

    $ python fig_appendix.py

Execution time: few hours
Author: Matheus Rolim Sales
Last modified: 19/02/2023
"""

import numpy as np # NumPy module
from joblib import Parallel, delayed # Module to create parallel loops
from functions import lyapunov, RTE_border, RTE_border_v2, RTE_border_v3 # Module with the standard map functions
from functions import corr_coef

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
# 
teps = np.linspace(1/100, 15/100, 29, endpoint=True)

# Obtains the largest Lyapunov exponent for the array above using parallel computing with Numba.
lyap = lyapunov(x, p, k, T)

# Obtains the rte for the array above using parallel computing.
# n_jobs=-1 uses all available threads. Set it to another value if you 
# wish to use less.
rte = [Parallel(n_jobs=-1)(delayed(RTE_border)(x, p, k[i], T, eps=teps[j]) for i in range(L)) for j in range(len(teps))]
rte2 = [Parallel(n_jobs=-1)(delayed(RTE_border_v2)(x, p, k[i], T, eps=teps[j]) for i in range(L)) for j in range(len(teps))]
rte3 = [Parallel(n_jobs=-1)(delayed(RTE_border_v3)(x, p, k[i], T, eps=teps[j]) for i in range(L)) for j in range(len(teps))]

with open("Data/cc_vs_eps.dat", "w") as df:
    cc1 = []
    cc2 = []
    cc3 = []
    for i in range(len(teps)):
        cc1.append(corr_coef(lyap, rte[i]))
        cc2.append(corr_coef(lyap, rte2[i]))
        cc3.append(corr_coef(lyap, rte3[i]))
        df.write("%.16f %.16f\n" % (teps[i], cc1[i], cc2[i], cc3[i]))