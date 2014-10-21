import numpy as np
from os import path

cdir = path.split(__file__)[0]
datadir = path.join(cdir, 'tables/')

eV = 1.60217657e-19  # [J]
c0 = 299792458  # [m/s]
h = 6.62606957e-34  # [m^2 kg / s]
kB = 1.3806488e-23  # [m^2 kg / s^2 / K]
T_CMB = 2.72548  # CMB temperature [K]


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
        self.data = {}  # dictionary (redshift : spectrum)

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical number density.
        """
        eps_data, n_data = self.data[z]
        return np.interp(eps, eps_data, n_data, 0, 0)

    def getEmin(self, z=0):
        """Minimum tabulated photon energy in [J]"""
        return self.data[z][0][0]

    def getEmax(self, z=0):
        """Maximum tabulated photon energy in [J]"""
        return self.data[z][0][-1]

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
        n = d[1] * 4 * np.pi / c0 / eps**2
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
            self.data[z] = eps, n

class EBL_Stecker05(EBL):
    name = 'IRB_Stecker05'
    info = 'cosmic infrared and optical background radiation model of Stecker at al. 2005'
    files = datadir+'IRO_Stecker05/z0'
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

class CRB_Biermann96(EBL):
    name = 'URB_Biermann96'
    info = 'cosmic radio background radiation model of Biermann et al. 1996'
    files = datadir+'CRB_Biermann96/z0'
    redshift = (0)

    def __init__(self):
        EBL.__init__(self)
        # d[0] : log(nu [GHz])
        # d[1-3] : log(I_nu [W / m^2 / Hz / sr])
        # radio galaxies, normal galaxies, normal galaxies (no evolution)
        d = np.genfromtxt(self.files, names=('nu', 'I1', 'I2', 'I3'))
        d.sort(order=['nu'], axis=0) # sort by frequency
        eps = 10**(d['nu']+9) * h
        n = (10**d['I1'] + 10**d['I2']) * 4*np.pi / c0 / h / eps
        self.data[0] = eps, n

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
    eps = logspace(-3, 1, 100) * eV
    ebl1 = EBL_Kneiske04()
    ebl2 = EBL_Kneiske10()
    ebl3 = EBL_Stecker05()
    ebl4 = EBL_Dole06()
    ebl5 = EBL_Franceschini08()

    figure()
    plot(eps/eV, ebl1.getDensity(eps), label="Kneiske '04")
    plot(eps/eV, ebl2.getDensity(eps), label="Kneiske '10 (lower limit)")
    plot(eps/eV, ebl3.getDensity(eps), label="Stecker '05")
    plot(eps/eV, ebl4.getDensity(eps), label="Dole '06")
    plot(eps/eV, ebl5.getDensity(eps), label="Franceschini '08")
    legend(loc='lower left')
    loglog()
    ylabel('$dn / d\epsilon$ [1/m$^3$/J]')
    xlabel('$\epsilon$ [eV]')
    savefig('IRO.png')
    show()


