"""Standard map functions

This module contains functions for the dynamical analysis of the standard map

    x_{n + 1} = x_{n} + p_{n + 1},
    p_{n + 1} = p_{n} - k * sin(x_{n}),

where n is the discrete time. It contains the following functions:

    * lyapunov - returns the largest Lyapunov exponent
    * stdmap - returns the time series 
    * RTE - returns the recurrence time entropy
    * FTRTE - returns the finite-time recurrence time entropy distribution
    * RTE_border - returns the recurrence time entropy considering border effects
    * FTRTE_border - returns the finite-time recurrence time entropy distribution considering border effects
    * get_trappingtimes - returns the trapping times
    * get_Qtau - returns the cumulative distribution of trapping times
    * plot_params - adjust the parameters for plotting and returns the color map used in Figs. 4 and 5

Author: Matheus Rolim Sales
Last modified: 02/12/2022
"""

import numpy as np # NumPy module
from numba import vectorize, njit # Numba module to create fast functions
from pyunicorn.timeseries import RecurrencePlot as RP # Pyunicorn module to create the recurrence plots
# Plotting module
import matplotlib.pyplot as plt
import matplotlib as mpl

@vectorize(['f8(f8, f8, f8, i8)', 'f4(f4, f4, f4, i4)'],
           target='parallel',
           nopython=True)
def lyapunov(x0, y0, k, N):
    """
    Calculate the largest Lyapunov exponent for the standard map given an initial condition (`x0`, `y0`).

    Parameters
    ------------
    x : float or array
        x-axis initial condition.
    y : float or array
        y-axis initial condition.
    k : float or array
        Non-linearity parameter of the standard map.
    N : int
        Number of iterations.

    Returns
    ------------
    out : float or array
        The largest Lyapunov exponent.
    """
    x = x0
    y = y0
    J = np.ones((2, 2), dtype=np.float64)
    beta0 = 0
    sumT11 = 0

    for i in range(N):
        # Iterates the map
        y = (y - k*np.sin(x)) % (2*np.pi)
        x = (x + y) % (2*np.pi)
        # Jacobian matrix
        J[0, 1] = -k*np.cos(x)
        J[1, 1] = 1.0 - k*np.cos(x)
        # Rotation angle
        beta = np.arctan((-J[1,0]*np.cos(beta0) + J[1,1]*np.sin(beta0))/(J[0,0]*np.cos(beta0) - J[0,1]*np.sin(beta0)))
        # Element 1,1 of the upper triangular matrix
        T11 = np.cos(beta0)*(J[0,0]*np.cos(beta) - J[1,0]*np.sin(beta)) - np.sin(beta0)*(J[0,1]*np.cos(beta) - J[1,1]*np.sin(beta))
        sumT11 += np.log(abs(T11))/np.log(2)
        # Update the rotation angle
        beta0 = beta
    
    lypnv = sumT11/N

    return lypnv
   
@njit
def stdmap(x0, y0, k, T):
    """
    Return the time series of the standard map given an initial condition (`x0`, `y0`).

    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    T : int
        The number of iterations (length of the orbit).

    Return
    ------
    out :  (T + 1, 2)-array
        The time series, where u[0:T + 1, 0] = x(t) and u[0:T + 1, 1] = y(t).
    """
    u = np.zeros((T + 1, 2))
    u[0, 0] = x0
    u[0, 1] = y0

    for i in range(1, T + 1):
        u[i, 1] = (u[i - 1, 1] - k*np.sin(u[i - 1, 0])) % (2*np.pi)
        u[i, 0] = (u[i - 1, 0] + u[i, 1]) % (2*np.pi)

        if u[i, 0] < 0: u[i, 0] += 2*np.pi
        if u[i, 0] > np.pi: u[i, 0] -= 2*np.pi
        if u[i, 1] < 0: u[i, 1] += 2*np.pi
        if u[i, 1] > np.pi: u[i, 1] -= 2*np.pi
        
    return u

def RTE(x0, y0, k, T, metric='supremum', return_last_pos=False):
    """
    Return the recurrence time entropy (rte) [1] given an initial condition (`x0`, `y0`).
    
    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    T : int
        The number of iterations (length of the orbit).
    metric : string (Optional)
        The metric for measuring distances in phase space. The valid options are "manhattan", "euclidean" or "supremum". (Default = "supremum").
    return_last_pos : bool (Optional)
        If True also return the last position of the orbit. (Default = False)

    Return
    ------
    out : tuple
        The white vertical entropy
        If return_last_pos = True, the second and third elements are last position in phase space (x, y).

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    """
    time_series = stdmap(x0, y0, k, T)
    rp = RP(time_series, metric=metric, normalize=False, threshold_std=10/100, silence_level=2)
    rte = rp.white_vert_entropy(w_min=1)
    if return_last_pos:
        return rte, time_series[:, 0][-1], time_series[:, 1][-1]
    else:
        return rte

def FTRTE(x0, y0, k, n, Ntot):
    """
    Return the distribuitions of the finite-time recurrence time entropy (ftrte) [1].

    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    n : int
        The number of iterations (length of the orbit - finite-time).
    Ntot : int
        The total number of iterations.

    Return
    ------
    out: tuple
        The distribuitions of the ftrte.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    """
    N = round(Ntot/n)
    ftrte = np.zeros(N)
    x = x0
    y = y0
    for i in range(N):
        ftrte[i], x, y = RTE(x, y, k, n, return_last_pos=True)

    return ftrte

@njit
def white_vertline_distr(recmat):
    """
    Returns the distribution of white vertical lines that not start nor end at the border of the RP.

    Parameters
    ----------
    recmat : 2D-array
        The recurrence matrix
    
    Returns
    -------
    out : 1D-array
        The distribution of white vertical lines
    """
    N = recmat.shape[0]
    P = np.zeros(N)
    for i in range(N):
        k = 0 # Counts the length of the white lines
        l = 0 # Checks if the white line is not on the border
        for j in range(N):
            if recmat[i, j] == 0 and l != 0:
                k += 1
            else:
                if k != 0:
                    P[k] += 1
                    k = 0
                elif recmat[i, j] == 1:
                    l = 1
    
    return P

def RTE_border(x0, y0, k, T, metric='supremum', lmin=1, return_last_pos=False):
    """
    Return the recurrence time entropy (rte) [1] given an initial condition (`x0`, `y0`) considering border effects when evaluating the distribution of white vertical lines [2].
    
    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    T : int
        The number of iterations (length of the orbit).
    metric : string (Optional)
        The metric for measuring distances in phase space. The valid options are "manhattan", "euclidean" or "supremum". (Default = "supremum").
    return_last_pos : bool (Optional)
        If True also return the last position of the orbit. (Default = False)

    Return
    ------
    out : tuple
        The white vertical entropy
        If return_last_pos = True, the second and third elements are last position in phase space (x, y).

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] Kraemer, K. H., Marwan, N. (2019): Border effect corrections for diagonal line based recurrence quantification analysis measures. - Physics Letters A, 383, 34, Art. 125977
    """
    time_series = stdmap(x0, y0, k, T)
    rp = RP(time_series, metric=metric, normalize=False, threshold_std=10/100, silence_level=2)
    recmat = rp.recurrence_matrix()
    p = white_vertline_distr(recmat)
    p = p[lmin:]
    p = np.extract(p != 0, p)

    p_normed = p/p.sum()

    if return_last_pos:
        return - (p_normed*np.log(p_normed)).sum(), time_series[-1, 0], time_series[-1, 1]
    else:
        return - (p_normed*np.log(p_normed)).sum()
    
def FTRTE_border(x0, y0, k, n, Ntot):
    """
    Return the distribuitions of the finite-time recurrence time entropy (ftrte) [1] considering border effects when evaluating the distribution of white vertical lines [2].

    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    n : int
        The number of iterations (length of the orbit - finite-time).
    Ntot : int
        The total number of iterations.

    Return
    ------
    out: tuple
        The distribuitions of the ftrte.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] Kraemer, K. H., Marwan, N. (2019): Border effect corrections for diagonal line based recurrence quantification analysis measures. - Physics Letters A, 383, 34, Art. 125977
    """
    N = round(Ntot/n)
    ftrte = np.zeros(N)
    x = x0
    y = y0
    for i in range(N):
        ftrte[i], x, y = RTE_border(x, y, k, n, return_last_pos=True)

    return ftrte

@njit
def get_trappingtimes(x, s0, sf):
    """
    Returns the "trapping times" given the RTE time series `x` and the trapping interval [`s0', `sf`].

    Parameters
    ----------

    x : 1D-array
        The RTE time series
    s0 : float
        Beginning of the interval
    sf : float
        Ending of the interval

    Returns
    -------
    out : 1D-list
        The trapping times
    """
    tau = []
    count = 0
    for i in range(len(x)):
        if x[i] >= s0 and x[i] <= sf:
            count += 1
        elif count != 0:
            tau.append(count)
            count = 0
    return tau

def get_Qtau(tau, nts):
    """
    Returns the cumulative distribuion of trapping times.

    Parameters
    ----------

    tau : 1D-array
        The trapping times
    nts : int
        Number of points of Q

    Returns
    -------
    
    out : (1D-array, 1D-array)
        The cumulative times and the cumulative distribuition of trapping times
    """
    t = np.logspace(np.log10(1), np.log10(max(tau)), nts)
    Q = np.zeros(len(t))
    N = len(tau)
    for j in range(len(t)):
        rec_tau = tau[np.where(tau > t[j])[0]]
        Ntau = len(rec_tau)
        Q[j] = Ntau/N

    return t, Q

def plot_params(fontsize=14, tick_labelsize=17, axes_labelsize=20, legend_fontsize=14):
    """
    Update the parameters of the plot.

    Return
    ------
    out : string
        The color map used in the colored plots.
    """
    plt.clf()
    plt.rc('font', size=fontsize)
    plt.rc('xtick', labelsize=tick_labelsize)
    plt.rc('ytick', labelsize=tick_labelsize)
    plt.rc('axes', labelsize=axes_labelsize)
    plt.rc('legend', fontsize=legend_fontsize)
    font = {'family' : 'stix'}
    plt.rc('font', **font)
    plt.rcParams["mathtext.fontset"] = "stix"
    cmap = 'nipy_spectral'
    mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally

    return cmap