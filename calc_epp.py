"""
Calculate the energy loss rate from electron pair production against isotropic photon fields
References:
(B70) Blumenthal 1970, Phys.Rev. D
(C92) Chodorowski et al. 1992, ApJ 400:181-185
"""

from numpy import *
from scipy import interpolate, integrate
import photonField


eV = 1.60217657e-19  # [J]
Mpc = 3.08567758e22  # [m]
c0 = 299792458.  # [m/s]
r0 = 2.817940e-15  # classical electron radius [m]
alpha = 7.297352e-3  # fine-structure constant
me = 9.10938291e-31  # electron mass [kg]
me_c2 = me * c0**2  # electron mass in [J/c^2]
mp = 1.67262178e-27  # proton mass [kg]


def lossRate(gamma, field, z=0):
    """
    Loss rate through electron pair production with the given photon background, cf. C92, equation 3.11
    gamma   : list of nucleus Lorentz factors
    field   : photon background
    z       : redshift
    Returns : 1/gamma dgamma/dx [1/Mpc]
    """

    _c = array([0.8048, 0.1459, 1.137e-3, -3.879e-6])
    _d = array([-86.07, 50.96, -14.45, 8/3.])
    _f = array([2.910, 78.35, 1837])
    def phi(k):
        """
        Parametrization of the integral 3.12 (C92)
        """
        # phi(k) for k < 25, eq. 3.14
        if k < 25:
            return pi/12 * (k-2)**4 / (1 + sum(_c * (k-2)**arange(1,5)))
        # phi(k) for k > 25, eq. 3.18 and 3.16
        phi = k * sum(_d * log(k)**arange(4))
        return phi / (1 - sum(_f * k**-arange(1,4)))

    def integrand(logk, gamma, field):
        """
        Integrand of equation 3.11 (C92), logarithmic version
        logk  : ln(k) = ln(2 gamma eps / (me c^2)), photon energy
        gamma : nucleus Lorentz factor
        field : photon background
        """
        k = exp(logk)
        eps = k * me_c2 / 2 / gamma  # photon energy [J] in lab frame
        n = field.getDensity(eps, z)  # spectral number density [1/m^3/J]
        n *= me_c2  # from substitution d eps / d k
        return n * phi(k) / k  # includes *k from substitution k -> ln(k)


    rate = zeros_like(gamma)
    err  = zeros_like(gamma)

    # minimum and maximum energy of the fields photons in units of me*c^2
    emin = field.getEmin() / me_c2
    emax = field.getEmax() / me_c2

    for i, g in enumerate(gamma):
        lkmin = log(max(2, 2*g*emin))
        lkmax = log(2*g * emax)
        lksep = linspace(lkmin, lkmax, 11)[1:-1]
        rate[i], err[i] = integrate.quad(integrand,
            lkmin, lkmax, points=lksep, args=(g, field))

    # prefactor of equation 3.11 (C92) and conversion [1/s] --> [1/Mpc]
    a = alpha * r0**2 * me / mp * Mpc
    return a * rate / gamma, a * err / gamma


# -------------------------------------------------
# Generate tables for the energy loss rate through pair production
# -------------------------------------------------
gamma = logspace(6, 14, 161)  # tabulated Lorentz factors
fname = 'data/epp_%s.txt'

fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12(),
    ]

for field in fields:
    print field.name
    rate = lossRate(gamma, field)[0]
    s = (rate > 1e-12)  # truncate if loss rate is < 10^-12 / Mpc

    fname  = 'data/epp_%s.txt' % field.name
    data   = c_[log10(gamma[s]), rate[s]]
    fmt    = '%.2f\t%.6e'
    header = 'Loss rate for electron-pair production with the %s\nlog10(gamma)\t1/gamma dgamma/dx [1/Mpc]' % field.info
    savetxt(fname, data, fmt=fmt, header=header)


# -------------------------------------------------
# Generate tables for the energy loss rate through pair production
# including time evolution of the photon field
# -------------------------------------------------
def gen_table_z(field, fname):
    redshifts = field.redshift  # list of available redshifts

    n = len(redshifts)
    m = len(gamma)
    rate = zeros((n, m))

    for i, z in enumerate(redshifts):
        print z
        rate[i] = lossRate(gamma, field, z)[0]

    # # truncate if loss rate is < 10^-12 / Mpc
    # s = (rate > 1e-12)
    table = (log10(gamma), rate)
    table = c_[gamma, rate.T]

    savetxt(fname, table, delimiter='\t',
        header='Loss rate for electron-pair production with the %s\nlog10(gamma)\t1/gamma dgamma/dx [1/Mpc]'%field.info)

# gamma = logspace(6, 14, 161)  # tabulated Lorentz factors
# fname = 'ElectronPairProduction/pair_rate_z_%s.txt'

# gen_table_z(photonField.CMB(), fname%'CMB')
# gen_table_z(photonField.KneiskeEBL(), fname%'KneiskeIRB')


# -------------------------------------------------
# Reformat CRPropa2 tables of differential spectrum of secondary electrons
# This should be reimplemented for extension to the other backgrounds,
# cross-checking and documentation.
# -------------------------------------------------
d1 = genfromtxt('tables/EPP/pair_spectrum_cmb.table', unpack=True)
d2 = genfromtxt('tables/EPP/pair_spectrum_cmbir.table', unpack=True)

# amplitudes dN/dEe(Ep)
A1 = d1[2].reshape((70, 170)) # CMB
A2 = d2[2].reshape((70, 170)) # CMB + IRB (which?)
A3 = A2 - A1  # IRB only

# normalize to 1
A1 = (A1.T / sum(A1, axis=1)).T
A3 = (A2.T / sum(A2, axis=1)).T

# save
savetxt('data/pair_spectrum_CMB.txt', A1, fmt='%.5e')
savetxt('data/pair_spectrum_IRB.txt', A3, fmt='%.5e')
