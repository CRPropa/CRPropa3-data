import numpy as np
from scipy import integrate


eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]


def calc_rate_eps(eps, xs, gamma, field, z=0):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    eps   : tabulated photon energies [J] in nucleus rest frame
    xs    : tabulated cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py

    Returns : inverse mean free path [1/Mpc]
    """
    F = integrate.cumtrapz(x=eps, y=eps*xs, initial=0)
    n = field.getDensity(np.outer(1./(2*gamma), eps), z)
    dx = np.mean(np.diff(np.log(eps)))  # value of log-spacing
    return integrate.romb(n * F / eps, dx=dx) / gamma * Mpc

def calc_rate_s(s_kin, xs, E, field):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    1/lambda = 1/(2E) * \int n((smax-m^2)/(4E)) / (smax-m^2) F(smax-m^2) dln(smax-m^2)
    F(smax-m^2) = \int_{smin-m^2}^{smax-m^2} sigma(s) (s - m^2) d(s-m^2)
    s = 2E eps (1 - cos(theta)) + m^2
    smax = 4E*eps + m^2
    smin = m^2

    s_kin : tabulated (s - m**2) for cross sections [J^2]
    xs    : tabulated cross sections [m^2]
    E     : (array of) cosmic ray energies [J]
    field : photon background, see photonField.py

    Returns : interaction rate [1/Mpc]
    """
    F = integrate.cumtrapz(x=s_kin, y=s_kin*xs, initial=0)
    n = field.getDensity(np.outer(1./(4*E), s_kin))
    ds = np.mean(np.diff(np.log(s_kin)))  # value of log-spacing
    return integrate.romb(n * F / s_kin, dx=ds) / 2 / E * Mpc

def calc_diffrate_eps(eps, xs, gamma, field):
    """
    Calculate the cumulative differential interaction rate against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    eps   : tabulated photon energies [J] in nucleus rest frame
    xs    : tabulated cross sections [m^2]
    E     : (array of) cosmic ray energies [J]
    field : photon background, see photonField.py

    Returns : cumulative differential interaction rates [1/Mpc]
    """
    # F = integrate.cumtrapz(x=s_kin, y=s_kin*xs, initial=0)
    # n = field.getDensity(np.outer(1./(4*E), s_kin))
    # y = n * F / s_kin**2 / E[:,np.newaxis] / 2 * Mpc
    # return integrate.cumtrapz(x=s_kin, y=y, initial=0)
    F = integrate.cumtrapz(x=eps, y=eps*xs, initial=0)
    n = field.getDensity(np.outer(1./(2*gamma), eps))
    y = n * F / eps**2 / gamma[:,np.newaxis] / *Mpc
    return integrate.cumtrapz(x=eps, y=y, initial=0)

def calc_diffrate_s(s_kin, xs, E, field):
    """
    Calculate the cumulative differential interaction rate against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    s_kin : tabulated (s - m**2) for cross sections [J^2]
    xs    : tabulated cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py

    Returns : cumulative differential interaction rates [1/Mpc]
    """
    F = integrate.cumtrapz(x=s_kin, y=s_kin*xs, initial=0)
    n = field.getDensity(np.outer(1./(4*E), s_kin))
    y = n * F / s_kin**2 / E[:,np.newaxis] / 2 * Mpc
    return integrate.cumtrapz(x=s_kin, y=y, initial=0)

def romb_truncate(x):
    """ Truncate array to largest size n = 2^i + 1 """
    i = int( np.floor(np.log2(n)) ) + 1
    return x[0:2**i+1]

def romb_pad_zero(x, n):
    """ Pad array with zeros """
    npad = n - len(x)
    return np.r_[x, np.zeros(npad)]

def romb_pad_logspaced(x, n):
    """ Pad array with log-linear increasing values """
    npad = n - len(x)
    dlx = np.mean( np.diff( np.log(x) ) )
    xpad = x[-1] * np.exp( dlx * np.arange(1, npad+1) )
    return np.r_[x, xpad]
