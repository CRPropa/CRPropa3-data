import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit


eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]

def invMFP(teps, txs, gamma, field):
    """
    Calculate the inverse mean free path for given tabulated cross sections against an isotropic photon background

    teps  : tabulated photon energies [J] in nucleus rest frame
    txs   : tabulated cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py

    Returns : inverse mean free path [1/Mpc]
    """
    # integral sigma(eps') eps' deps'  in nucleus rest frame
    tF1 = integrate.cumtrapz(x=teps, y=txs*teps, initial=0)
    # substitute k = ln(eps/J)  photon energy in observer frame
    kmin, kmax = np.log((field.getEmin(), field.getEmax()))
    tk = np.linspace(kmin, kmax, 5000)  # upsample
    exptk = np.exp(tk)

    rate = np.zeros_like(gamma)

    for i, g in enumerate(gamma):
        # integrand: n(eps) * F1(2 * g * eps) / eps
        tf2 = field.getDensity(exptk) / exptk * np.interp(2*g*exptk, teps, tF1)

        # limit integration to range where integrand != 0
        innz = tf2.nonzero()[0]
        if len(innz) < 2:
            continue  # rate = 0
        k0, k1 = tk[innz[0]], tk[innz[-1]]

        ksep = np.linspace(k0, k1, 11)[1:-1]  # subdivisions

        rate[i] = integrate.quad(
            lambda k: np.interp(k, tk, tf2),
            k0, k1, points=ksep)[0]

    return rate * Mpc / (2 * gamma**2)


def invMFP_fast(eps, xs, gamma, field):
    """
    Calculate the inverse mean free path for given tabulated
    cross sections against an isotropic photon background:
    Fully vectorized version using Romberg integration.

    The tabulated cross sections need to be of length n = 2^i + 1
    and the tabulation points log-linearly spaced
    eps   : photon energies [J] in nucleus rest frame
    xs    : cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py

    Returns : inverse mean free path [1/Mpc]
    """
    F = integrate.cumtrapz(x=eps, y=eps*xs, initial=0) / eps
    n = field.getDensity( np.outer(1 / (2 * gamma), eps) )
    dx = np.mean(np.diff(np.log(eps)))  # average log-spacing
    return integrate.romb(n * F, dx=dx) * Mpc / gamma


def romb_truncate(x):
    """ Truncate to largest size n = 2^i + 1 """
    i = int( np.floor(np.log2(n)) ) + 1
    return x[0:2**i+1]

def romb_pad_zero(x, n):
    """ Pad with zeros """
    npad = n - len(x)
    return np.r_[x, np.zeros(npad)]

def romb_pad_logspaced(x, n):
    """ Pad with log-linear increasing values """
    npad = n - len(x)
    dlx = np.mean( np.diff( np.log(x) ) )
    xpad = x[-1] * np.exp( dlx * np.arange(1, npad+1) )
    return np.r_[x, xpad]

def romb_pad_expo(x, y, x2, nfit=10):
    """ Pad with exponential extrapolation using last nfit points """
    # check if extrapolation is possible
    if y[-1] == 0:
        return  np.r_[y, np.zeros(len(x2)-len(x))]

    # y = a * x^b  --> log(y) = log(a) + b * log(x)
    lx = np.log(x[-nfit:])
    ly = np.log(y[-nfit:])

    def f(x, la, b):
        return la + b * x

    la, b = curve_fit(f, lx, ly, p0 = (0, ly[0]))[0]

    y2 = np.exp( la + b * np.log(x2) )
    y2[0:len(y)] = y
    return y2
