import numpy as np
from scipy.integrate import cumulative_trapezoid, romb, quad
import os

eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]


def calc_rate_eps(eps, xs, gamma, field, z=0, cdf=False):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    eps   : tabulated photon energies [J] in nucleus rest frame
    xs    : tabulated cross sections [m^2]
    gamma : (array of) nucleus Lorentz factors
    field : photon background, see photonField.py
    z     : redshift
    cdf   : calculate cumulative differential rate

    Returns :
        interaction rate 1/lambda(gamma) [1/Mpc] or
        cumulative differential rate d(1/lambda)/d(s_kin) [1/Mpc/J^2]
    """
    F = cumulative_trapezoid(x=eps, y=eps * xs, initial=0)
    n = field.getDensity(np.outer(1. / (2 * gamma), eps), z)
    if cdf:
        y = n * F / eps**2
        return cumulative_trapezoid(x=eps, y=y, initial=0) / np.expand_dims(gamma, -1) * Mpc
    else:
        y = n * F / eps
        dx = mean_log_spacing(eps)
        return romb(y, dx=dx) / gamma * Mpc


def calc_rate_s(s_kin, xs, E, field, z=0, cdf=False):
    """
    Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
    The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

    s_kin : tabulated (s - m**2) for cross sections [J^2]
    xs    : tabulated cross sections [m^2]
    E     : (array of) cosmic ray energies [J]
    field : photon background, see photonField.py
    z     : redshift
    cdf   : calculate cumulative differential rate

    Returns :
        interaction rate 1/lambda(gamma) [1/Mpc] or
        cumulative differential rate d(1/lambda)/d(s_kin) [1/Mpc/J^2]
    """

    if cdf:
        # precalculate the field integral if it not exisist and load it afterwards
        calculateDensityIntegral(field)
        file = "data/fieldDensity/" + field.name + ".txt"
        densityIntegral = np.loadtxt(file)

        # interpolate
        I = np.zeros((len(E), len(s_kin)))
        for j in range(len(E)):
            I[j,:] = np.interp(s_kin/ 4 / E[j], densityIntegral[:,0], densityIntegral[:,1])

        # calculate cdf
        y = np.array([xs * s_kin for i in range(len(E))]) * I
        return cumulative_trapezoid(y = y, x = s_kin, initial=0) / 8 / np.expand_dims(E, -1)**2 * Mpc    
    else:
        F = cumulative_trapezoid(x=s_kin, y=s_kin * xs, initial=0)
        n = field.getDensity(np.outer(1. / (4 * E), s_kin), z)
        y = n * F / s_kin
        ds = mean_log_spacing(s_kin)
        return romb(y, dx=ds) / 2 / E * Mpc

def calculateDensityIntegral(field):
    """ Precalculate the integral over the density 
        int_{Emin}^{Emax} n(eps) / eps^2  deps 
        and save as a file.

        field : photon background, see photonField.py
    """

    # precalc the photon density integral 
    Emax = field.getEmax()
    Emin =  1e4 / 4 / 1e23 * eV # min(s_kin) / 4 / max(E_e)
    alpha = np.logspace(np.log10(Emin), np.log10(Emax), 10000) # lower boundary of the integral.

    # check if file already exist
    folder = "data/fieldDensity/"
    if not os.path.isdir(folder):
        os.makedirs(folder)
    file = folder + field.name + ".txt"
    if os.path.isfile(file):
        return # file already existst no calculation necessary

    # calculate integral
    I_gamma = np.zeros_like(alpha)
    for i in range(len(alpha)):
        I_gamma[i] = quad(lambda E: field.getDensity(E) / E**2, a = alpha[i], b = Emax, full_output=1)[0]

    # save file
    header = "# calculate integral n(e)/e^2 de from eMin to eMax, where eMax is the maximal photon energy of the background \n"
    try: 
        git_hash = gh.get_git_revision_hash()
        header += "Produced with crpropa-data version: "+git_hash+"\n"
        header += "# eMin [eV]\tintegral\n"
    except:
        header += "# eMin [eV]\tintegral\n"
    data = np.c_[alpha, I_gamma]
    fmt = '%.4e\t%8.7e'
    np.savetxt(file, data, fmt = fmt, header = header)


def mean_log_spacing(x):
    """ <Delta log(x)> """
    return np.mean(np.diff(np.log(x)))


def romb_truncate(x, n):
    """ Truncate array to largest size n = 2^i + 1 """
    i = int(np.floor(np.log2(n))) + 1
    return x[0:2**i + 1]


def romb_pad_zero(x, n):
    """ Pad array with zeros """
    npad = n - len(x)
    return np.r_[x, np.zeros(npad)]


def romb_pad_logspaced(x, n):
    """ Pad array with log-linear increasing values """
    npad = n - len(x)
    dlx = np.mean(np.diff(np.log(x)))
    xpad = x[-1] * np.exp(dlx * np.arange(1, npad + 1))
    return np.r_[x, xpad]
