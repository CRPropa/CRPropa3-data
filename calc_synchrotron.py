import os
import numpy as np
from scipy.special import kv
from scipy import integrate
import gitHelp as gh


def synchrotron_spectrum(xval):
    """
    Calculate cumulative synchrotron spectrum.
    Follows:
      J.~D. Jackson, Classical Electrondynamics.
      Wiley, 3rd ed., p. 681, eq. 14.91

    F(x) = (\int_{x}^{\infinity} K_{5/3}(x') dx') 
    for x = E_gamma / E_critical, where E_critical = hbar * 3/2 * c/rho * (E/(mc**2))**2
    E_gamma : Energy of synchrotron photon
    E       : Energy of particle
    rho     : gyroradius 

    Returns : 
      The cumulative synchrotron function
    """
    F = np.zeros(len(xval))
    for i, value in enumerate(xval):
        a = xval[i:]
        F[i] = integrate.trapz(x = a, y = kv(5. / 3., a))
    for i,value in enumerate(xval):
        b = integrate.cumtrapz(x = xval, y = xval * F, initial = 0)
    return b / b[-1]


def compute_spectrum(x, outputName):
    """
    Cumulative differential synchrotron spectrum. 
    This implementation follows:
      J.~D. Jackson, Classical Electrondynamics.
      Wiley, 3rd ed., p. 681, eq. 14.91

    Input
    . x: fraction between photon frequency and critical frequency
    . outputName: name of output file
    """
    cdf = synchrotron_spectrum(x)
    lx = np.log10(x)
    data = np.c_[lx, cdf]
    # Add git hash of crpropa-data repository to header
    try:
        git_hash = gh.get_git_revision_hash()
        header = 'Produced with crpropa-data version: '+git_hash+'\nx\t: photon frequency to critical frequency fraction\nlog10(x)\tCDF\n'
    except: 
        header = 'x\t: photon frequency to critical frequency fraction\nlog10(x)\tCDF\n'
    fmt = '%3.2f\t%7.6e'
    np.savetxt(outputName, data, fmt = fmt, header = header)

def plot(specFile, plotFile):
    """
    Make simple plot for sanity checks.

    Input
    . specFile: file containing the synchrotron spectrum
    """
    data = np.loadtxt(specFile)
    x = 10 ** data[:, 0]
    y = data[:, 1]

    y = np.diff(y) / np.diff(x)
    x = 10 ** ((np.log10(x[:-1]) + np.diff(np.log10(x))) / 2.)
    # y = np.diff(y) / np.diff(data[])

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(x, y)
    plt.loglog()
    plt.grid()
    plt.xlim(1e-7, 1e5)
    plt.ylim(1e-6, 10.)
    plt.xlabel('$x \\equiv \\nu / \\nu_c$')
    plt.ylabel('$f(x)$')
    plt.savefig(plotFile)


# ----------------------------------------------------------------
# main
# ----------------------------------------------------------------
if __name__ == '__main__':
    
    x = np.logspace(-10., 4, 1401)

    plotDir = 'plots'
    resDir = 'data/Synchrotron'
    if not os.path.exists(plotDir):
            os.makedirs(plotDir)
    if not os.path.exists(resDir):
        os.makedirs(resDir)
    outputName  = '%s/spectrum.txt' % resDir
    
    compute_spectrum(x, outputName)

    plotName = '%s/sync.png' % plotDir
    plot(outputName, plotName)


