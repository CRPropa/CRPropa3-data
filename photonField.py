import numpy as np
from os import path

cdir = path.split(__file__)[0]
datadir = path.join(cdir, 'tables/')

eV = 1.60217657e-19  # [J]
erg = 1e-7  # [J]
c0 = 299792458  # [m/s]
h = 6.62606957e-34  # [m^2 kg / s]
kB = 1.3806488e-23  # [m^2 kg / s^2 / K]
T_CMB = 2.72548  # CMB temperature [K]

# --------------------------------------------------------
# interfaces
# --------------------------------------------------------
class CMB:
    """
    Cosmic microwave background radiation
    """
    name = 'CMB'
    info = 'CMB'
    redshift = None
    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical number density.
        """
        return 8*np.pi / c0**3 / h**3 * eps**2 / (np.exp(eps/(kB*T_CMB)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV

class EBL:
    """
    Base class for extragalactic background light (EBL) models
    """
    def __init__(self):
        self.data = {}  # dictionary {redshift : (eps, dn/deps)}

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical spectral number density.
        The tabulated data is interpolated linearly in log-log; no extrapolation is performed.
        """
        tab_eps, tab_n = self.data[z]
        # return np.interp(eps, tab_eps, tab_n, 0, 0)  # linear interpolation
        x  = np.log10(eps)
        tx = np.log10(tab_eps)
        ty = np.log10(tab_n)
        return 10**np.interp(x, tx, ty, -np.inf, -np.inf)  # log-log interpolation

    def getEmin(self, z=0):
        """Minimum tabulated photon energy in [J]"""
        return self.data[z][0][0]

    def getEmax(self, z=0):
        """Maximum tabulated photon energy in [J]"""
        return self.data[z][0][-1]

# --------------------------------------------------------
# EBL (optical and infrared) models
# --------------------------------------------------------
class EBL_Kneiske04(EBL):
    name = 'IRB_Kneiske04'
    info = 'cosmic infrared and optical background radiation model of Kneiske et al. 2004'
    files = datadir+'IRO_Kneiske04/all_z'
    redshift = np.linspace(0, 5, 51)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : eps [eV]
        # d[1-51] : n(eps), [1/m^3/eV]
        d = np.genfromtxt(self.files, unpack=True)
        eps = d[0] * eV
        n = d[1:] / eV
        for i,z in enumerate(self.redshift):
            self.data[z] = eps, n[i]

class EBL_Kneiske04_old(EBL):
    name = 'IRB_Kneiske04'
    info = 'cosmic infrared and optical background radiation model of Kneiske et al. 2004'
    files = datadir+'IRO_Kneiske04/old/%.1f'
    redshift = (0, 0.2, 0.4, 0.6, 1, 2, 3, 4)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # x : wavelength in [mu m]
            # y : lambda I_lambda [nW/m^2/sr]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            wl = (10**x * 1e-6)  # wavelength in [m]
            eps = h * c0 / wl
            n = (10**y * 1e-9) * 4*np.pi/c0 / eps**2
            self.data[z] = eps[::-1], n[::-1]

class EBL_Kneiske10(EBL):
    name = 'IRB_Kneiske10'
    info = 'cosmic infrared and optical background radiation lower limit model of Kneiske et al. 2010'
    files = datadir+'IRO_Kneiske10/%.1f'
    redshift = (0, .1, .3, .8, 2)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # x : wavelength in [mu m]
            # y : lambda I_lambda [nW/m^2/sr]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            wl = (10**x * 1e-6)  # wavelength in [m]
            eps = h * c0 / wl
            n = (10**y * 1e-9) * 4*np.pi/c0 / eps**2
            self.data[z] = eps[::-1], n[::-1]

class EBL_Dole06(EBL):
    name = 'IRB_Dole06'
    info = 'cosmic infrared and optical background radiation model of Dole et al. 2006'
    files = datadir+'IRO_Dole06/0.0'
    redshift = (0)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : lambda [mu m]
        # d[1] : n(eps), [W/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # photon energy [J]
        n = d[1] * 4*np.pi / c0 / eps**2
        self.data[0] = eps[::-1], n[::-1]

class EBL_Franceschini08(EBL):
    name = 'IRB_Franceschini08'
    info = 'cosmic infrared and optical background radiation model of Franceschini et al. 2008'
    files = datadir+'IRO_Franceschini08/%1.1f'
    redshift = (0, .2, .4, .6, .8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # x : log10(eps / eV)
            # y : eps * dn/deps [1/cm^3]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            eps = 10**x * eV
            n = 10**y / eps * 1e6
            n /= (1 + z)**3  # make comoving
            self.data[z] = eps, n

class EBL_Stecker05_old(EBL):
    # data file obtained from O. Kalashevs propation code, "does not include UV part"
    name = 'IRB_Stecker05'
    info = 'cosmic infrared and optical background radiation model of Stecker at al. 2005'
    files = datadir+'IRO_Stecker05/data1.txt'
    redshift = np.linspace(0, 5, 26)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : log(eps/Hz)
        # d[1-26] : eps*n(eps), not comoving
        d = np.genfromtxt(self.files, unpack=True)
        eps = 10**d[0] * h
        n = d[1:] / eps
        for i,z in enumerate(self.redshift):
            # convert n(eps) to comoving density
            self.data[z] = eps, n[i] / (1+z)**3

class EBL_Stecker05(EBL):
    name = 'IRB_Stecker05'
    info = 'cosmic infrared and optical background radiation model of Stecker at al. 2005'
    files = datadir+'IRO_Stecker05/data2.txt'
    redshift = np.linspace(0, 5, 26)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : log10(eps/eV)
        # d[1-26] : log10(eps*n(eps)/cm^3), not comoving
        d = np.genfromtxt(self.files, unpack=True)
        eps = 10**d[0] * eV
        n = 10**d[1:] / eps * 1e6
        for i,z in enumerate(self.redshift):
            # convert n(eps) to comoving density
            self.data[z] = eps, n[i] / (1+z)**3

class EBL_Finke10(EBL):
    name = 'IRB_Finke10'
    info = 'cosmic infrared and optical background radiation model of Finke et al. 2010 (Model C)'
    files = datadir+'IRO_Finke10/z%.2f.dat'
    redshift = np.arange(0, 5, 0.01)

    def __init__(self):
        EBL.__init__(self)
        for z in self.redshift:
            # d[0] : eps / eV
            # d[1] : comoving energy density in erg / cm^3
            d = np.genfromtxt(self.files % z, unpack=True)
            eps = d[0] * eV
            n   = d[1] * erg * 1e6 / eps**2 # [J/m^3]
            self.data[z] = eps, n

class EBL_Gilmore12(EBL):
    name = 'IRB_Gilmore12'
    info = 'cosmic infrared and optical background radiation model of Gilmore et al. 2012 (Evolving dust model, arXiv:1104.0671)'
    files = datadir+'IRO_Gilmore12/eblflux_fiducial.dat'
    redshift = np.array([0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])

    def __init__(self):
        EBL.__init__(self)
        # d[0] : rest frame wavelength in [angstrom]
        # d[1-21] : proper flux [erg/s/cm^2/ang/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-10)  # [J]
        n = d[1:] * 1e-7 / c0 / 1e-4 * d[0] * 4*np.pi / eps**2
        for i,z in enumerate(self.redshift):
            n[i] /= (1 + z)**3  # make comoving
            self.data[z] = eps[::-1], n[i][::-1]

class EBL_Dominguez11(EBL):
    name = 'IRB_Dominguez11'
    info = 'cosmic infrared and optical background radiation model of Dominguez et al. 2011 (arXiv:1007.1459)'
    files = datadir+'IRO_Dominguez11/ebl_dominguez11.out'
    redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])

    def __init__(self):
        EBL.__init__(self)
        # d[0] : rest frame wavelength in [mu m]
        # d[1-18] : proper flux [nW/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # [J]
        n = d[1:] * 1e-9 / c0 * 4*np.pi / eps**2
        for i,z in enumerate(self.redshift):
            self.data[z] = eps[::-1], n[i][::-1]  # sort by ascending energy

class EBL_Dominguez11_upper(EBL):
    name = 'IRB_Dominguez11_upper'
    info = 'cosmic infrared and optical background radiation model of Dominguez et al. 2011 (arXiv:1007.1459)'
    files = datadir+'IRO_Dominguez11/ebl_upper_uncertainties_dominguez11.out'
    redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])

    def __init__(self):
        EBL.__init__(self)
        # d[0] : rest frame wavelength in [mu m]
        # d[1-18] : proper flux [nW/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # [J]
        n = d[1:] * 1e-9 / c0 * 4*np.pi / eps**2
        for i,z in enumerate(self.redshift):
            self.data[z] = eps[::-1], n[i][::-1]  # sort by ascending energy

class EBL_Dominguez11_lower(EBL):
    name = 'IRB_Dominguez11_lower'
    info = 'cosmic infrared and optical background radiation model of Dominguez et al. 2011 (arXiv:1007.1459)'
    files = datadir+'IRO_Dominguez11/ebl_lower_uncertainties_dominguez11.out'
    redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])

    def __init__(self):
        EBL.__init__(self)
        # d[0] : rest frame wavelength in [mu m]
        # d[1-18] : proper flux [nW/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h * c0 / (d[0] * 1e-6)  # [J]
        n = d[1:] * 1e-9 / c0 * 4*np.pi / eps**2
        for i,z in enumerate(self.redshift):
            self.data[z] = eps[::-1], n[i][::-1]  # sort by ascending energy

# --------------------------------------------------------
# CRB (radio) models
# --------------------------------------------------------
class CRB_Protheroe96:
    """
    Universal Radio Background from Protheroe & Bierman 1996.
    Taken from EleCa: Parametrization from fit to plots from paper.
    """
    name = "CRB_Protheroe96"
    info = 'cosmic radio background model of Protheroe & Biermann 1996'

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        if z != 0:
            print 'CRB_Protheroe96: no evolution available!'

        p0 = -2.23791e+01
        p1 = -2.59696e-01
        p2 =  3.51067e-01
        p3 = -6.80104e-02
        p4 =  5.82003e-01
        xmin = -6.003761  # 4.1e-12 eV
        xmax = -0.315515  # 2e-6 eV

        x = np.log10(np.r_[eps] / h / 1e9)
        I = p0 + p1 * x + p2 * x**2 + p3 * x**3 / (np.exp(p4 * x) - 1)
        I = 4 * np.pi / (h * c0) * (10**I / eps)
        I[x < xmin] = 0
        I[x > xmax] = 0
        return I

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 4.1e-12 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 2e-6 * eV

class CRB_ARCADE2:
    def getDensity(self, eps):
        """
        Spectral number density dn/deps [1/m^3/J] at z = 0.
        """
        # T = 1.26 +- 0.09 K (nu/GHz)^-2.6 +- 0.04, see Holder 2012
        T = 1.26 * (nu/1e9)**-2.6
        return 8*np.pi / c0**3 / h**3 * eps**2 / (np.exp(eps/(kB*T)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV


if __name__ == '__main__':
    from pylab import *
    eps = logspace(-3, 1, 200) * eV
    x  = eps / eV
    y1 = EBL_Kneiske04().getDensity(eps) * eps
    y3 = EBL_Stecker05().getDensity(eps) * eps
    y5 = EBL_Franceschini08().getDensity(eps) * eps
    y6 = EBL_Finke10().getDensity(eps) * eps
    y7 = EBL_Dominguez11().getDensity(eps) * eps
    y8 = EBL_Gilmore12().getDensity(eps) * eps
    upper = EBL_Dominguez11_upper().getDensity(eps) * eps
    lower = EBL_Dominguez11_lower().getDensity(eps) * eps

    figure()
    plot(x, y1, label="Kneiske '04")
    plot(x, y3, label="Stecker '05")
    plot(x, y5, label="Franceschini '08")
    plot(x, y6, label="Finke '10")
    plot(x, y7, 'k', label="Dominguez '11")
    fill_between(x, lower, upper, edgecolor='none', facecolor='grey')
    plot(x, y8, label="Gilmore '12")
    legend(loc='lower left')
    loglog()
    ylim(1e1, 1e6)
    ylabel('$\epsilon ~ dn/d\epsilon$ [1/m$^3$]')
    xlabel('$\epsilon$ [eV]')
    savefig('figures/IRO.png')
