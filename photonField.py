import numpy as np
from crpropa import eV, erg, c_light, h_planck, k_boltzmann, hertz, ccm
import os
import gitHelp as gh
import pandas as pd

cdir = os.path.split(__file__)[0]
datadir = os.path.join(cdir, 'tables/')

class PhotonField(object):
    """Base class for photon fields"""

    def __init__(self):
        self.name = 'PhotonField'
        self.info = 'Base class photon field'
        self.energy = [] #[eV]
        self.redshift = None
        self.photonDensity = [] #[eV^-1 cm^-3]
        self.outdir = 'data/Scaling'

    def createFiles(self):
        try:
            git_hash = gh.get_git_revision_hash()
            addHash = True
        except:
            addHash = False

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        with open(self.outdir + "/" + self.name + "_photonEnergy.txt", 'w') as f:
            f.write('# '+self.info+'\n')
            if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
            f.write("# photon energies in [J]\n")
            for e in self.energy:
                f.write("{}\n".format(e * eV))  # [J]
        if self.redshift is not None:
            with open(self.outdir + "/" + self.name + "_redshift.txt", 'w') as f:
                f.write('# '+self.info+'\n')
                if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
                f.write("# redshift\n")
                for z in self.redshift:
                    f.write("{}\n".format(np.round(z, 2)))
        with open(self.outdir + "/" + self.name + "_photonDensity.txt", 'w') as f:
            f.write('# '+self.info+'\n')
            if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
            f.write("# Comoving photon number density in [m^-3], format: d(e1,z1), ... , d(e1,zm), d(e2,z1), ... , d(e2,zm), ... , d(en,zm)\n")
            for i, densSlice in enumerate(self.photonDensity):
                #Including redshift evolution
                try:
                    for d in densSlice:
                        f.write("{}\n".format(d * self.energy[i] / ccm))  # [# / m^3], comoving
                #When no redshift is included the densSlice is a 1d array
                except TypeError:
                    f.write("{}\n".format(densSlice * self.energy[i] / ccm))  # [# / m^3], comoving
        print("done: " + self.name)

# --------------------------------------------------------
# interfaces
# --------------------------------------------------------
class CMB(PhotonField):
    """
    Cosmic microwave background radiation
    """
    
    def __init__(self):
        super(CMB, self).__init__()
        self.name = 'CMB'
        self.info = 'Cosmic Microwave Background, T_CMB = 2.72548 K'
        self.T_CMB = 2.72548  # CMB temperature [K]
        self.energy = np.logspace(-10, -1, 101) # [eV]
        self.photonDensity = self.getDensity(self.energy * eV) * (eV * ccm) # [1/eVcm^3]

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical number density.
        """
        return 8*np.pi / c_light**3 / h_planck**3 * eps**2 / np.expm1(eps / (k_boltzmann * self.T_CMB)) 

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV


class EBL(PhotonField):
    """
    Base class for extragalactic background light (EBL) models
    """
    def __init__(self):
        super(EBL, self).__init__()
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

    def getEnergy(self, z=0):
        return self.data[z][0][:]

# --------------------------------------------------------
# EBL (optical and infrared) models
# --------------------------------------------------------
class EBL_Kneiske04(EBL):
    """ IRB model from Kneiske 2004 """

    def __init__(self):
        super(EBL_Kneiske04, self).__init__()
        self.name = 'IRB_Kneiske04'
        self.info = 'cosmic infrared and optical background radiation model of Kneiske et al. 2004'
        self.files = datadir + 'EBL_Kneiske_2004/all_z'
        self.redshift = np.linspace(0, 5, 51)

        # d[0] : eps [eV]
        # d[1-51] : n(eps), [1/m^3/eV]
        d = np.genfromtxt(self.files, unpack=True)
        eps = d[0] * eV
        n = d[1:] / eV
        for i,z in enumerate(self.redshift):
            self.data[z] = eps, n[i]

        d = np.genfromtxt(self.files)
        self.photonDensity = []
        self.energy = []
        for i,fieldSlice in enumerate(d):
            for j,entry in enumerate(fieldSlice):
                if j != 0:
                    d[i][j] /= 1e6  # 1/eVm^3 -> 1/eVcm^3
            self.photonDensity.append(fieldSlice[1:])
            self.energy.append(fieldSlice[0])
        

class EBL_Kneiske10(EBL):
    """ IRB model from Kneiske 2010 """

    def __init__(self):
        super(EBL_Kneiske10, self).__init__()
        self.name = 'IRB_Kneiske10'
        self.info = 'cosmic infrared and optical background radiation lower limit model of Kneiske et al. 2010'
        self.files = datadir + 'EBL_Kneiske_2010/%.1f'
        self.redshift = (0, .1, .3, .8, 2)
        for z in self.redshift:
            # x : wavelength in [mu m]
            # y : lambda I_lambda [nW/m^2/sr]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            wl = (10**x * 1e-6)  # wavelength in [m]
            eps = h_planck * c_light / wl
            n = (10**y * 1e-9) * (4 * np.pi / c_light) / eps**2
            self.data[z] = eps[::-1], n[::-1]
    
    #TODO: Check if necessary information is available to create photonDensity and energy arrays
    #NOTE: EBL_Kneiske10 was and is still *not* available in CRPropa. If this should be changed also 
    #      CRPropa's PhotonBackground.h needs to be updated.
    def createFiles(self):
        print("WARNING "+self.__class__.__name__+": Not all relevant data available. No Files produced")
        return

class EBL_Dole06(EBL):
    """ IRB model from Dole 2006"""

    def __init__(self):
        super(EBL_Dole06, self).__init__()
        self.name = 'IRB_Dole06'
        self.info = 'cosmic infrared and optical background radiation model of Dole et al. 2006'
        self.files = datadir + 'EBL_Dole_2006/0.0'
        self.redshift = (0)
        # d[0] : lambda [mu m]
        # d[1] : n(eps), [W/m^2/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h_planck * c_light / (d[0] * 1e-6)  # photon energy [J]
        n = d[1] * (4 * np.pi / c_light) / eps**2
        self.data[0] = eps[::-1], n[::-1]

    #TODO: Check if necessary information is available to create photonDensity and energy arrays
    #NOTE: EBL_Dole06 was and is still *not* available in CRPropa. If this should be changed also 
    #      CRPropa's PhotonBackground.h needs to be updated.
    def createFiles(self):
        print("WARNING "+self.__class__.__name__+": Not all relevant data available. No Files produced")
        return

class EBL_Franceschini08(EBL):
    """ IRB model from Fanceschini 2008 """

    def __init__(self):
        super(EBL_Franceschini08, self).__init__()
        self.name = 'IRB_Franceschini08'
        self.info = 'cosmic infrared and optical background radiation model of Franceschini et al. 2008'
        self.files = datadir + 'EBL_Franceschini_2008/%1.1f'
        self.fileDir = datadir + 'EBL_Franceschini_2008/'
        self.redshift = np.linspace(0., 2., 11)
        for z in self.redshift:
            # x : log10(eps / eV)
            # y : eps * dn/deps [1/cm^3]
            x, y = np.genfromtxt(self.files%z, unpack=True)
            eps = 10**x * eV
            n = 10**y / eps * 1e6
            n /= (1 + z)**3  # make comoving
            self.data[z] = eps, n
        self.loadData()
        self.loadEnergy()
        self.loadPhotonDensity()


    def loadData(self):
        """Load data from all redshift files into a single dictionary"""
        fileList = os.listdir(self.fileDir)
        self.fieldData = {}
        # read in data. Note that there is a different energy range for each redshift!
        for i, file in enumerate(sorted(fileList)):
            if "README.txt" not in file:
                data = np.genfromtxt(self.fileDir + file, unpack=True)
                eps = data[0]  # [log10eV]
                dens = data[1]  # [log10(eps*dn/deps in cm^-3)]
                self.fieldData[np.round(self.redshift[i], 2)] = (eps, dens)

    def loadEnergy(self):
        """determine global min and max eps and then generate a logspace from them."""
        epsMin, epsMax = [], []
        for key in self.fieldData.keys():
            epsMin.append(min(self.fieldData[key][0]))
            epsMax.append(max(self.fieldData[key][0]))
        self.energy = 10 ** np.linspace(min(epsMin), max(epsMax), 41)  # [eV]

    def loadPhotonDensity(self):
        """interpolate all redshift-dependent photon field densities on this range"""
        photonField = []
        for key in self.fieldData.keys():
            dens = np.interp(self.energy, 10 ** self.fieldData[key][0], self.fieldData[key][1])
            dens = 10 ** dens / self.energy  # log10(e*dn/de in 1/cm^3) -> 1/(eVcm^3)
            photonField.append(dens)
        photonField = np.swapaxes(np.array(photonField), 0, 1)
        self.photonDensity = [zSlice / (1 + self.redshift)**3 for zSlice in photonField]  # make comoving

class EBL_Stecker05(EBL):
    """ IRB model from Stecker 2005 """

    def __init__(self):
        super(EBL_Stecker05, self).__init__()
        self.name = 'IRB_Stecker05'
        self.info = 'cosmic infrared and optical background radiation model of Stecker at al. 2005'
        self.files = datadir + 'EBL_Stecker_2005/data2.txt'
        self.redshift = np.linspace(0, 5, 26)
        # d[0]    : log10(eps/eV)
        # d[1-26] : log10(eps*n(eps)/cm^3), not comoving
        d = np.genfromtxt(self.files, unpack=True)
        eps = 10**d[0] * eV
        n = 10**d[1:] / eps * 1e6
        for i,z in enumerate(self.redshift):
            # convert n(eps) to comoving density
            self.data[z] = eps, n[i] / (1+z)**3

        data = d.transpose()
        self.energy = []
        self.photonDensity = []
        for i, zSlice in enumerate(data):
            eps = 10**zSlice[0]  # [eV]
            self.energy.append(eps)
            dens = 10**zSlice[1:] / eps  # [1/eVcm^3]
            dens /= (self.redshift + 1)**3   # make comoving
            self.photonDensity.append(dens)

class EBL_Finke10(EBL):
    """ IRB model from Finke 2010"""

    def __init__(self):
        super(EBL_Finke10, self).__init__()
        self.name = 'IRB_Finke10'
        self.info = 'cosmic infrared and optical background radiation model of Finke et al. 2010 (Model C)'
        self.files = datadir + 'EBL_Finke_2010/z%.2f.dat'
        self.redshift = np.arange(0, 5, 0.01)
        for z in self.redshift:
            # d[0] : eps / eV
            # d[1] : comoving energy density in erg / cm^3
            d = np.genfromtxt(self.files % z, unpack=True)
            eps = d[0] * eV
            n   = d[1] * erg * 1e6 / eps**2 # [J/m^3]
            self.data[z] = eps, n

        self.fileDir = datadir + "EBL_Finke_2010/"
        self.fileList = os.listdir(self.fileDir)
        d = pd.DataFrame()
        col = 0
        for file in sorted(self.fileList):
            if "README.txt" not in file:
                data = np.genfromtxt(self.fileDir + file, unpack=True)
                eps = data[0]  # [eV] It's the same energy for all redshifts.
                data = pd.DataFrame(data[1],columns=[col/100])
                col += 1
                d = pd.concat([d,data], join='outer', axis=1)
        self.photonDensity = []
        self.energy = []
        for i,e in enumerate(eps):
            #TODO: Check if the factor 6.2415091e11 is correct
            #      1/eV * centimeter**3 = 6.2415091e12
            #      Where does the additional factor 0.1 come from?
            dens = np.array(list(d.iloc[i])) * 6.2415091e11 / e**2  # [1/eVcm^3]
            self.photonDensity.append(dens)
            self.energy.append(e)

class EBL_Gilmore12(EBL):
    """ IRB model from Gilmore 2012 """

    def __init__(self):
        super(EBL_Gilmore12, self).__init__()
        self.name = 'IRB_Gilmore12'
        self.info = 'cosmic infrared and optical background radiation model of Gilmore et al. 2012 (Evolving dust model, arXiv:1104.0671)'
        self.files = datadir + 'EBL_Gilmore_2012/eblflux_fiducial.dat'
        self.redshift = np.array([0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6,\
                                 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])
        # d[0] : rest frame wavelength in [angstrom]
        # d[1-21] : proper flux [erg/s/cm^2/ang/sr]
        d = np.genfromtxt(self.files, unpack=True)
        eps = h_planck * c_light / (d[0] * 1e-10)  # [J]
        n = d[1:] * erg / 1e-4 * d[0] / eps**2 * (4 * np.pi / c_light) 
        for i,z in enumerate(self.redshift):
            n[i] /= (1 + z)**3  # make comoving
            self.data[z] = eps[::-1], n[i][::-1]

        wavelength = d[0]  # angstrom
        eps = 12.39842e3 / wavelength  # [eV]
        self.photonDensity = []
        self.energy = []
        d = pd.DataFrame(d)
        for i in range(len(eps)):
            fieldSlice = np.array(list(d[i])[1:]) * 4 * np.pi / (100 * c_light) * wavelength[i] * erg  # eV/cm^3
            fieldSlice /= eps[i]**2  # 1/eVcm^3
            self.photonDensity.append(fieldSlice / eV)  # /eV?!
            self.energy.append(eps[i])
        # invert, because lambda is antiprop to energy
        self.photonDensity = [x / (1 + self.redshift)**3 for x in reversed(self.photonDensity)]  # also make comoving
        self.energy = [e for e in reversed(self.energy)]

class EBL_Dominguez11(EBL):
    """ IRB model from Dominguez 2011 
    
    EBL intensities from the paper "Extragalactic background light inferred from AEGIS galaxy-SED-type fractions", A. Dominguez et al., 2011, MNRAS, 410, 2556
    """

    def __init__(self, which='best'):
        """ Constructor
        
        Input:
          which : \"best\" for the best fit model, \"upper\" or \"lower\" for the upper/lower uncertaincy.
        """
        super(EBL_Dominguez11, self).__init__()
    
        if which == 'best':
            fname = 'EBL_Dominguez_2011/ebl_dominguez11.out'
            self.name = 'IRB_Dominguez11'
        elif which == 'upper':
            fname = 'EBL_Dominguez_2011/ebl_upper_uncertainties_dominguez11.out'
            self.name = 'IRB_Dominguez11_upper'
        elif which == 'lower':
            fname = 'EBL_Dominguez_2011/ebl_lower_uncertainties_dominguez11.out'
            self.name = 'IRB_Dominguez11_lower'
        else:
            raise ValueError('EBL_Dominguez11 only provides "best", "upper" and "lower" models')

        self.info = 'cosmic infrared and optical background radiation ({}) model of Dominguez et al. 2011 (arXiv:1007.1459)'.format(which)
        self.redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])

        # d[0] : rest frame wavelength in [mu m]
        # d[1-18] : proper flux [nW/m^2/sr]
        d = np.genfromtxt(datadir + fname, unpack=True)
        eps = h_planck * c_light / (d[0] * 1e-6)  # [J]
        n = d[1:] * 1e-9 / eps**2 * (4 * np.pi / c_light)
        for i, z in enumerate(self.redshift):
            self.data[z] = eps[::-1], n[i][::-1]  # sort by ascending energy

        d = np.genfromtxt(datadir + fname, unpack=False)
        eps = [h_planck*c_light / (eV * ccm) / fieldSlice[0]*eV for fieldSlice in d]  # [eV] | 1.238842 = h*c/1µ eV
        #                          nW->W : J->eV : 1/m^3->1/cm^3 : 1/sm^2sr->1/m^3
        n = np.array([fieldSlice[1:] * 1e-9 * eV / 1e6 / eps[i]**2 * (4 * np.pi / c_light) for i, fieldSlice in enumerate(d)])  # [1/eVcm^3] 
        energy = []
        for i,x in enumerate(n):
            energy.append(eps[i]/eV)
        self.photonDensity = [x for x in reversed(n)]
        self.energy = [e for e in reversed(energy)]

class EBL_Stecker16(EBL):
    """ IRB model from Stecker 2016 """

    def __init__(self, which='upper'):
        """
        Constructor 
        
        Input
          which : \"upper\" or \"lower\" for an estimation of the uncertaincy
        """
        super(EBL_Stecker16, self).__init__()

        if which == 'upper':
            fname = 'EBL_Stecker_2016/comoving_enerdens_up.csv'
        elif which == 'lower':
            fname = 'EBL_Stecker_2016/comoving_enerdens_lo.csv'
        else:
            raise ValueError('EBL_Stecker16 only provides "upper" and "lower" models')

        self.name = 'IRB_Stecker16_%s' % which
        self.info = 'cosmic infrared and optical background radiation (%s) model of Stecker et al. 2016 (arXiv:1605.01382)' % which
        self.redshift = np.linspace(0, 5, 501)
        
        # rows:    log10(photon energy / eV) ranging from -2.84 to 1.14 in steps of 0.01
        # columns: redshifts from 0 to 5 in steps of 0.01
        # units:   erg / Hz / cm^3 --> 
        d = np.genfromtxt(datadir + fname, delimiter=',').T
        eps = 10**np.arange(-2.84, 1.14001, 0.01) * eV
        n = d * erg / h_planck * 1E6 / eps
        for i, z in enumerate(self.redshift):
            self.data[z] = eps, n[i]

        d_t = d.T
        eps_eV = eps / eV # [eV]
        nu = eps / c_light  # [Hz]
        self.photonDensity = []
        self.energy = []
        for i, dens in enumerate(d_t):
            dens = dens * 6.2415091e11 * nu[i] / eps_eV[i]**2  # 1/eVcm^3
            self.photonDensity.append(dens)
            self.energy.append(eps_eV[i])

# --------------------------------------------------------
# CRB (radio) models
# --------------------------------------------------------
class URB_Protheroe96(PhotonField):
    """
    Universal Radio Background from Protheroe & Bierman 1996.
    Taken from EleCa implementation.

    Reference
    R. J. Protheroe and P. L. Biermann
    Astroparticle Physics 6 (1996) 45.
    """

    def __init__(self):

        super(URB_Protheroe96, self).__init__()
        self.name = "URB_Protheroe96"
        self.info = "Universal Radio Background from Protheroe & Bierman Astroparticle Physics 6 (1996) 45."
        self.redshift = None

        logEmin = np.log10(self.getEmin())
        logEmax = np.log10(self.getEmax())
        eps = np.logspace(logEmin, logEmax, 101)
        self.photonDensity = np.array([self.getDensity(e) for e in eps]) / (ccm**-1  * eV**-1)
        self.energy = eps / eV

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        p0 = -2.23791e+01
        p1 = -2.59696e-01
        p2 = 3.51067e-01
        p3 = -6.80104e-02
        p4 = 5.82003e-01
        p5 = -2.00075e+00
        p6 = -1.35259e+00
        p7 = -7.12112e-01  # xbreak

        eps = np.r_[eps]
        x = np.log10(eps / h_planck / 1e9)
        I = p0 + p1 * x + p2 * x**2 + p3 * x**3 / (np.exp(p4 * x) - 1)
        I[x > p7] += p6 + p5 * x[x > p7] - p2 * x[x > p7]**2
        I = 4 * np.pi / (h_planck * c_light) * (10**I / eps)

        I[eps < self.getEmin()] = 0
        I[eps > self.getEmax()] = 0
        return I

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 4.1e-12 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 2E-6 * eV # 0.825e-6 * eV

class URB_Fixsen11(PhotonField):
    """
    Universal Radio Background as measured by ARCADE2.
    Note that the frequency range in this reference is more narrow than for other models.
    Therefore, this should be used carefully.

    Reference:
      D. J. Fixsen et al.
      The Astrophysical Journal 734 (2011) 5.
      https://arxiv.org/abs/0901.0555
    """

    def __init__(self):
        super(URB_Fixsen11, self).__init__()
        self.name = 'URB_Fixsen11'
        self.info = 'Universal Radio Background as measured by ARCADE2 D. J. Fixsen et al. ApJ 734 (2011) 5'
        self.redshift = None
        self.T_CMB = 2.72548

        logEmin = np.log10(self.getEmin())
        logEmax = np.log10(self.getEmax())
        eps = np.logspace(logEmin, logEmax, 101)
        self.photonDensity = np.array([self.getDensity(e) for e in eps]) / (ccm**-1 * eV**-1)
        self.energy = eps / eV

    def getDensity(self, eps, z = 0.):
        """
        Spectral number density dn/deps [1/m^3/J] at z = 0.
        """
        eps = np.r_[eps]
        nu = eps / h_planck
        T = self.T_CMB + 24.1 * np.power(nu / 3.1e8, -2.6)
        I = 8. * np.pi / c_light ** 3 / h_planck ** 3 * eps ** 2 / (np.expm1(eps / (k_boltzmann * T)))
        I[eps < self.getEmin()] = 0.
        I[eps > self.getEmax()] = 0.
        return I

    def getEmin(self, z = 0.):
        """Minimum effective photon energy in [J]"""
        return 2.2e6 * hertz * h_planck

    def getEmax(self, z = 0.):
        """Maximum effective photon energy in [J]"""
        return 1e10 * hertz * h_planck

class URB_Nitu21(PhotonField):
    """
    Universal Radio Background from Nitu et al. 2021.
    Reference:
      I. C. Nitu, H. T. J. Bevings, J. D. Bray, A. M. M. Scaife
      Astroparticle Physics 126 (2021) 102532.
      https://arxiv.org/abs/2004.13596
    """

    def __init__(self):
        super(URB_Nitu21, self).__init__()
        self.name = 'URB_Nitu21'
        self.info = 'Universal Radio Background from Nitu et al. Astropart. Phys. 126 (2021) 102532'
        self.redshift = None

        logEmin = np.log10(self.getEmin())
        logEmax = np.log10(self.getEmax())
        eps = np.logspace(logEmin, logEmax, 101)
        self.photonDensity = np.array([self.getDensity(e) for e in eps]) / (ccm**-1  * eV**-1)
        self.energy = eps / eV

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        p0 = -1.9847e1
        p1 = -2.9857e-1
        p2 = -2.6984e-1
        p3 = 9.5394e-2
        p4 = -4.9059e-2
        p5 = 4.4297e-3
        p6 = 7.6038e-3
        p7 = -1.9690e-3
        p8 = -2.2573e-4
        p9 = 1.1762e-4
        p10 = -9.9443e-6
        p = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]

        eps = np.r_[eps]
        nu = eps / h_planck
        I = 0.
        for k in range(len(p)):
            I += (p[k] * np.power(np.log10(nu / 1e6), k))
        I = 10. ** I
        I = 4 * np.pi / (h_planck * c_light) * (I / eps)

        I[eps < self.getEmin()] = 0.
        I[eps > self.getEmax()] = 0.

        return I

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e3 * hertz * h_planck

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 1e12 * hertz * h_planck

# --------------------------------------------------------
#   main:   plot comparison of IRB models
# --------------------------------------------------------
if __name__ == '__main__':
    from pylab import *
    eps = logspace(-3, 1, 200) * eV
    x  = eps / eV
    c =  eps**2 / eV
    y1   = c * EBL_Kneiske04().getDensity(eps)
    y3   = c * EBL_Stecker05().getDensity(eps)
    y5   = c * EBL_Franceschini08().getDensity(eps)
    y6   = c * EBL_Finke10().getDensity(eps)
    y7   = c * EBL_Dominguez11().getDensity(eps)
    y8   = c * EBL_Gilmore12().getDensity(eps)
    y7up = c * EBL_Dominguez11('upper').getDensity(eps)
    y7lo = c * EBL_Dominguez11('lower').getDensity(eps)
    y9up = c * EBL_Stecker16('upper').getDensity(eps)
    y9lo = c * EBL_Stecker16('lower').getDensity(eps)

    figure()
    plot(x, y1, label='Kneiske 2004')
    plot(x, y3, label='Stecker 2005')
    plot(x, y5, label='Franceschini 2008')
    plot(x, y6, label='Finke 2010')
    plot(x, y7, label='Dominguez 2011')
    plot(x, y8, label='Gilmore 2012')
    fill_between(x, y7lo, y7up, facecolor='m', edgecolor='none', alpha=0.2, zorder=-1, label='Dominguez 2011 (limits)')
    fill_between(x, y9lo, y9up, facecolor='g', edgecolor='none', alpha=0.2, zorder=-1, label='Stecker 2016 (limits)')

    legend(loc='lower center', fontsize='x-small')
    loglog()
    grid()
    ylabel('$\epsilon^2 ~ dn/d\epsilon$ [eV/m$^3$]')
    xlabel('$\epsilon$ [eV]')
    savefig('plots/EBL.png')
    show()
