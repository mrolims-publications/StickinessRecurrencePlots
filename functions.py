"""Standard map functions

This module contains functions for the dynamical analysis of the standard map

    x_{n + 1} = x_{n} + p_{n + 1},
    p_{n + 1} = p_{n} - k * sin(x_{n}),

where n is the discrete time, used in the publication: "Stickiness and recurrence plots: an entropy based approach".

It contains the following functions:

    * lyapunov - return the largest Lyapunov exponent
    * stdmap - return the time series 
    * RTE - return the recurrence time entropy
    * FTRTE - return the finite-time recurrence time entropy distribution
    * RTE_border - return the recurrence time entropy considering border effects
    * FTRTE_border - return the finite-time recurrence time entropy distribution considering border effects
    * get_trappingtimes - return the trapping times
    * get_Qtau - return the cumulative distribution of trapping times
    * plot_params - adjust the parameters for plotting and returns the color map used in Figs. 4 and 5

Author: Matheus Rolim Sales
Last modified: 08/02/2023
"""

import numpy as np # NumPy module
from numba import vectorize, njit # Numba module to create fast functions
from pyunicorn.timeseries import RecurrencePlot as RP # Pyunicorn module to create the recurrence plots
# Plotting modules
import matplotlib.pyplot as plt
import matplotlib as mpl

@vectorize(['f8(f8, f8, f8, i8)', 'f4(f4, f4, f4, i4)'],
           target='parallel',
           nopython=True)
def lyapunov(x0, y0, k, N):
    """
    Calculate the Lyapunov exponent of the standard map map.

    This function uses a vectorized implementation to optimize performance.

    Parameters
    ----------
    x0 : float
        The initial x-coordinate of the map.
    y0 : float
        The initial y-coordinate of the map.
    k : float
        The constant parameter of the map.
    N : int
        The number of iterations of the map.

    Returns
    -------
    lyapunov_exponent : float
        The Lyapunov exponent of the map.

    Examples
    --------
    >>> x0 = 0.1
    >>> y0 = 0.1
    >>> k = 0.1
    >>> N = 10000
    >>> lyapunov(x0, y0, k, N)
    0.1...

    >>> x0 = [0.1, 0.2, 0.3]
    >>> y0 = 0.1
    >>> k = 0.1
    >>> N = 10000
    >>> lyapunov(x0, y0, k, N)
    [0.1..., 0.2..., 0.3...]

    >>> x0 = [0.1, 0.2, 0.3]
    >>> y0 = 0.1
    >>> k = [0.1, 0.5, 1.0]
    >>> N = 10000
    >>> lyapunov(x0, y0, k, N)
    [[0.1..., 0.2..., 0.3...], [0.4..., 0.5..., 0.6...], [0.7..., 0.8..., 0.9...]]

    Notes
    -----
    The Lyapunov exponent is a measure of the exponential rate of separation of nearby trajectories in a dynamical system. In this function, the standard map is iterated `N` times, and the Lyapunov exponent is calculated as the average of the logarithm of the magnitude of the upper triangular matrix (Jacobian matrix) element 1,1 divided by `N` (Eckmann-Ruelle method [1]).

    References
    ----------

    [1] J.-P. Eckmann and D. Ruelle, “Ergodic theory of chaos and strange attractors,” Reviews Modern Physics 57, 617–656 (1985).

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
def stdmap(x0, y0, k, N):
    """
    Iterate the standard map.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    x0 : float
        The initial x-coordinate of the map.
    y0 : float
        The initial y-coordinate of the map.
    k : float
        The constant parameter of the map.
    N : int
        The number of iterations of the map.

    Returns
    -------
    map_iterates : ndarray, shape (N+1, 2)
        The iterates of the standard map, including the initial state.

    Examples
    --------
    >>> x0 = 0.1
    >>> y0 = 0.1
    >>> k = 0.1
    >>> N = 10000
    >>> u = stdmap(x0, y0, k, N)

    Notes
    -----
    The standard map is a two-dimensional, nonlinear, area-preserving map that is often used as a simple model for Hamiltonian systems. In this function, the map is iterated `N` times, with the result stored in the `map_iterates` array. The values of x and y are limited to the interval [0, 2pi].

    """
    u = np.zeros((N + 1, 2))
    u[0, 0] = x0
    u[0, 1] = y0

    for i in range(1, N + 1):
        u[i, 1] = (u[i - 1, 1] - k*np.sin(u[i - 1, 0])) % (2*np.pi)
        u[i, 0] = (u[i - 1, 0] + u[i, 1]) % (2*np.pi)

        if u[i, 0] < 0: u[i, 0] += 2*np.pi
        if u[i, 0] > np.pi: u[i, 0] -= 2*np.pi
        if u[i, 1] < 0: u[i, 1] += 2*np.pi
        if u[i, 1] > np.pi: u[i, 1] -= 2*np.pi
        
    return u

def RTE(x0, y0, k, T, metric='supremum', return_last_pos=False):
    """
    Return the recurrence time entropy (RTE) [2, 3] for the standard map using pyunicorn [1].

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
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').
    return_last_pos : bool, optional
        Whether to return the final x, y position of the orbit. The default is False.
    
    Returns
    -------

    out : float or tuple of float, float, float
        If return_last_pos is False, returns the RTE.
        If return_last_pos is True, returns a tuple of the RTE, final x-coordinate and final y-coordinate.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    """
    time_series = stdmap(x0, y0, k, T)
    rp = RP(time_series, metric=metric, normalize=False, threshold_std=10/100, silence_level=2)
    rte = rp.white_vert_entropy(w_min=1)
    if return_last_pos:
        return (rte, time_series[:, 0][-1], time_series[:, 1][-1])
    else:
        return rte

def FTRTE(x0, y0, k, n, Ntot, metric='supremum'):
    """
    Return the finite-time recurrence time entropy (FTRTE) [1-3] distibution.

    Parameters
    ----------
    x0 : float
        The initial value of the x-coordinate.
    y0 : float
        The initial value of the y-coordinate.
    k : float
        The non-linearity parameter of the map.
    n : int
        The time window (finite-time).
    Ntot : int
        The total number of iterations.
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').

    Return
    ------
    out: narray
        Array with the FTRTE distribution.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    """
    N = round(Ntot/n)
    ftrte = np.zeros(N)
    x = x0
    y = y0
    for i in range(N):
        ftrte[i], x, y = RTE(x, y, k, n, metric=metric, return_last_pos=True)

    return ftrte

@njit
def white_vertline_distr(recmat):
    """
    Calculate the distribution of the lengths of white vertical lines in a binary matrix.

    Parameters
    ----------
    recmat : numpy.ndarray
        A 2-dimensional binary numpy array (recurrence matrix).

    Returns
    -------
    numpy.ndarray
        An array containing the count of white vertical lines for each length.
        
    Examples
    --------
    >>> recmat = np.array([[0, 0, 0, 1, 1], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0]])
    >>> white_vertline_distr(recmat)
    array([3., 2., 1.])

    Notes
    -----
    The input binary matrix is assumed to be a square matrix of size N x N and the lines on the border are excluded [1].

    References
    ----------
    [1] K. H. Kraemer and N. Marwan, Border effect corrections for diagonal line based recurrence quantification analysis measures, Physics Letters A 383, 125977 (2019)
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

def RTE_border(x0, y0, k, T, metric='supremum', lmin=1, eps=10/100, return_last_pos=False):
    """
    Return the recurrence time entropy (RTE) [1-3] given an initial condition (x0, y0) considering border effects [4] when evaluating the distribution of white vertical lines.

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
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').
    lmin : int, optional
        Minimal length of white vertical lines used in the RTE computation (default=1).
    teps : float, optional
        Threshold for the recurrence plot in unit of time series standard deviation (default=10/100).
    return_last_pos : bool, optional
        If True, also return the last position of the orbit. (default=False).

    Returns
    -------
    out : tuple or float
        If `return_last_pos=False` (default), returns the RTE. 
        If `return_last_pos=True`, returns a tuple (RTE, x, y), where `x` and `y` are the last position in phase space.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    [4] K. H. Kraemer and N. Marwan, Border effect corrections for diagonal line based recurrence quantification analysis measures, Physics Letters A 383, 125977 (2019)
    """
    time_series = stdmap(x0, y0, k, T)
    rp = RP(time_series, metric=metric, normalize=False, threshold_std=eps, silence_level=2)
    recmat = rp.recurrence_matrix()
    p = white_vertline_distr(recmat)
    p = p[lmin:]
    p = np.extract(p != 0, p)

    p_normed = p/p.sum()

    if return_last_pos:
        return (- (p_normed*np.log(p_normed)).sum(), time_series[-1, 0], time_series[-1, 1])
    else:
        return - (p_normed*np.log(p_normed)).sum()
    
def RTE_border_v2(x0, y0, k, T, metric='supremum', lmin=1, eps=10/100, return_last_pos=False):
    """
    Return the recurrence time entropy (RTE) [1-3] given an initial condition (x0, y0) considering border effects [4] when evaluating the distribution of white vertical lines.

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
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').
    lmin : int, optional
        Minimal length of white vertical lines used in the RTE computation (default=1).
    teps : float, optional
        Threshold for the recurrence plot in unit of time series standard deviation (default=10/100).
    return_last_pos : bool, optional
        If True, also return the last position of the orbit. (default=False).

    Returns
    -------
    out : tuple or float
        If `return_last_pos=False` (default), returns the RTE. 
        If `return_last_pos=True`, returns a tuple (RTE, x, y), where `x` and `y` are the last position in phase space.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    [4] K. H. Kraemer and N. Marwan, Border effect corrections for diagonal line based recurrence quantification analysis measures, Physics Letters A 383, 125977 (2019)
    """
    time_series = stdmap(x0, y0, k, T)
    eps = max(np.std(time_series[:, 0]), np.std(time_series[:, 1]))*eps
    rp = RP(time_series, metric=metric, normalize=False, threshold=eps, silence_level=2)
    recmat = rp.recurrence_matrix()
    p = white_vertline_distr(recmat)
    p = p[lmin:]
    p = np.extract(p != 0, p)

    p_normed = p/p.sum()

    if return_last_pos:
        return (- (p_normed*np.log(p_normed)).sum(), time_series[-1, 0], time_series[-1, 1])
    else:
        return - (p_normed*np.log(p_normed)).sum()
    
def RTE_border_v3(x0, y0, k, T, metric='supremum', lmin=1, eps=10/100, return_last_pos=False):
    """
    Return the recurrence time entropy (RTE) [1-3] given an initial condition (x0, y0) considering border effects [4] when evaluating the distribution of white vertical lines.

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
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').
    lmin : int, optional
        Minimal length of white vertical lines used in the RTE computation (default=1).
    teps : float, optional
        Threshold for the recurrence plot in unit of time series standard deviation (default=10/100).
    return_last_pos : bool, optional
        If True, also return the last position of the orbit. (default=False).

    Returns
    -------
    out : tuple or float
        If `return_last_pos=False` (default), returns the RTE. 
        If `return_last_pos=True`, returns a tuple (RTE, x, y), where `x` and `y` are the last position in phase space.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    [4] K. H. Kraemer and N. Marwan, Border effect corrections for diagonal line based recurrence quantification analysis measures, Physics Letters A 383, 125977 (2019)
    """
    time_series = stdmap(x0, y0, k, T)
    eps = np.sqrt(np.std(time_series[:, 0])**2 + np.std(time_series[:, 1])**2)*eps
    rp = RP(time_series, metric=metric, normalize=False, threshold=eps, silence_level=2)
    recmat = rp.recurrence_matrix()
    p = white_vertline_distr(recmat)
    p = p[lmin:]
    p = np.extract(p != 0, p)

    p_normed = p/p.sum()

    if return_last_pos:
        return (- (p_normed*np.log(p_normed)).sum(), time_series[-1, 0], time_series[-1, 1])
    else:
        return - (p_normed*np.log(p_normed)).sum()
    
def FTRTE_border(x0, y0, k, n, Ntot, metric='supremum', lmin=1, teps=10/100):
    """
    Calculates the finite-time recurrence time entropy (FTRTE) [1-3] considering border effects [4] on the distribution of white vertical lines.

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
    metric : str, optional
        The metric for measuring distances in phase space. Possible values are 'supremum', 'manhattan', 'euclidean' (default='supremum').
    lmin : int, optional
        Minimal length of white vertical lines used in the RTE computation (default=1).
    teps : float, optional
        Threshold for the recurrence plot in unit of time series standard deviation (default=10/100).

    Returns
    -------
    out : ndarray
        Array with the FTRTE distribution.
    """

    """
    Return the finite-time recurrence time entropy (FTRTE) [1-3] distibution considering border effects when evaluating the distribution of white vertical lines [4].

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
    out: narray
        Array with the FTRTE distribution.

    References
    ----------
    [1] http://www.pik-potsdam.de/~donges/pyunicorn/api/timeseries/recurrence_plot.html
    [2] M. A. Little et al., Exploiting nonlinear recurrence and fractal scaling properties for voice disorder detection, BioMedical Engineering OnLine 6, 23 (2007)
    [3] K. H. Kraemer et al., Recurrence threshold selection for obtaining robust recurrence characteristics in different embedding dimensions, Chaos 28, 085720 (2018)
    [4] K. H. Kraemer and N. Marwan, Border effect corrections for diagonal line based recurrence quantification analysis measures, Physics Letters A 383, 125977 (2019)
    """
    N = round(Ntot/n)
    ftrte = np.zeros(N)
    x = x0
    y = y0
    for i in range(N):
        ftrte[i], x, y = RTE_border(x, y, k, n, metric=metric, lmin=lmin, teps=teps, return_last_pos=True)

    return ftrte

@njit
def get_trappingtimes(x, s0, sf):
    """
    Calculates the trapping times in the RTE time series `x` based on the trapping interval [`s0`, `sf`].

    Parameters
    ----------
    x : ndarray
        1D numpy array representing the RTE time series.
    s0 : float
        The start of the trapping interval.
    sf : float
        The end of the trapping interval.

    Returns
    -------
    out : ndarray
        1D numpy array of trapping times.
    """
    tau = []
    count = 0
    for i in range(len(x)):
        if x[i] >= s0 and x[i] <= sf:
            count += 1
        elif count != 0:
            tau.append(count)
            count = 0
    return np.array(tau)

def get_Qtau(tau, nts):
    """
    Calculates the cumulative distribution of trapping times.

    Parameters
    ----------
    tau : ndarray
        An array of trapping times.
    nts : int
        The number of points to use for the cumulative distribution.

    Returns
    -------
    t : ndarray
        An array of the cumulative times.
    Q : ndarray
        An array of the cumulative distribution of trapping times.
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

    Returns
    -------
    cmap : string
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

def corr_coef(x, y):
    """
    Calculate the Pearson correlation coefficient between two arrays x and y.

    Parameters
    ----------
    x : ndarray
        The first input array. Must have the same number of elements as `y`.
    y : ndarray
        The second input array. Must have the same number of elements as `x`.

    Returns
    -------
    correlation_coefficient : float
        The Pearson correlation coefficient between `x` and `y`.

    Examples
    --------
    >>> x = [1, 2, 3, 4, 5]
    >>> y = [2, 4, 6, 8, 10]
    >>> corr_coef(x, y)
    1.0

    >>> x = [1, 2, 3, 4, 5]
    >>> y = [10, 8, 6, 4, 2]
    >>> corr_coef(x, y)
    -1.0

    >>> x = [1, 2, 3, 4, 5]
    >>> y = [2, 4, 7, 8, 10]
    >>> corr_coef(x, y)
    0.957...

    Notes
    -----
    The correlation coefficient is a value between -1 and 1, indicating the strength and direction of the linear relationship between `x` and `y`. A value of -1 indicates a perfect negative linear relationship, a value of 1 indicates a perfect positive linear relationship, and a value of 0 indicates no linear relationship.

    """
    std_x = np.std(x)
    std_y = np.std(y)
    covxy = np.cov(x, y)[0][1]
    cc = covxy/(std_x*std_y)

    return cc