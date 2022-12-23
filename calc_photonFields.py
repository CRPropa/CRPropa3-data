"""
    purpose: convert all native CRPropa photon field data from the
    CRPropa-data repo to a unified format for CRPropa to read.

    The data format for photon fields should be the same for each field:
    - a photonEnergy file with photon energies [e1, e2, ...en] in J
    - if there is redshift-dependence: a redshift file [z1, z2, zm]
    - a photonDensity file with the comoving number density with format
      [d(e1,z1), ... , d(e1,zm), d(e2,z1), ... , d(e2,zm), ... , d(en,zm)] in 1/m^3

    usage: this file should be placed in the CRPropa-data/ folder which you cloned from the CRPropa-data repo
    output: 3 data files per photon field, encoding photonEnergy, redshift, photonDensity
"""

import numpy as np
import pandas as pd 
import os
import gitHelp as gh
import photonField as pf
import crpropa as crp

cm3 = crp.centimeter**3 # [m^3]

def IRB_Stecker05(fileDir, outDir):
    name = 'IRB_Stecker05'
    info = '# cosmic infrared and optical background radiation model of Stecker at al. 2005'
    redshift = np.linspace(0., 5., 26)
    filePath = fileDir + "EBL_Stecker_2005/data2.txt"
    data = np.genfromtxt(filePath)
    energy = []
    photonField = []
    for i, zSlice in enumerate(data):
        eps = 10**zSlice[0]  # [eV]
        energy.append(eps)
        dens = 10**zSlice[1:] / eps  # [1/eVcm^3]
        dens /= (redshift + 1)**3   # make comoving
        photonField.append(dens)
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Gilmore12(fileDir, outDir):
    name = "IRB_Gilmore12"
    info = "# These tables contain the data for the background flux and associated optical depths of gamma rays for the WMAP5+Fixed ('fixed') and Evolving Dust ('fiducial') models presented in Gilmore, Somerville, Primack, and Dominguez (2012), ArXiv:1104.0671v2"
    filePath = fileDir + "EBL_Gilmore_2012/eblflux_fiducial.dat"
    redshift = np.array([0.0, 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])
    d = np.genfromtxt(filePath, unpack=True)
    wavelength = d[0]  # angstrom
    eps = 12.39842e3 / wavelength  # [eV]
    photonField = []
    energy = []
    d = pd.DataFrame(d)
    for i in range(len(eps)):
        fieldSlice = np.array(list(d[i])[1:]) * 4 * np.pi / (100 * crp.c_light) * wavelength[i] * crp.erg  # eV/cm^3
        fieldSlice /= eps[i]**2  # 1/eVcm^3
        photonField.append(fieldSlice / crp.eV)  # /eV?!
        energy.append(eps[i])
    # invert, because lambda is antiprop to energy
    photonField = [x / (1 + redshift)**3 for x in reversed(photonField)]  # also make comoving
    energy = [e for e in reversed(energy)]
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Finke10(fileDir, outDir):
    name = "IRB_Finke10"
    redshift = np.round(np.linspace(0., 4.99, 500), 2)
    info = "# Extragalactic background light model from Finke et al. 2010, DOI:10.1088/0004-637X/712/1/238, Files obtained from http://www.phy.ohiou.edu/~finke/EBL/"
    fileDir = fileDir + "EBL_Finke_2010/"
    fileList = os.listdir(fileDir)
    d = pd.DataFrame()
    col = 0
    for file in sorted(fileList):
        if "README.txt" not in file:
            data = np.genfromtxt(fileDir + file, unpack=True)
            eps = data[0]  # [eV]
            data = pd.DataFrame(data[1],columns=[col/100])
            index = data.index.copy()
            col += 1
            #d = pd.concat([d,data], axis=1, join_axes=[data.index])
            d = pd.concat([d,data], join='outer', axis=1)
            #print(d.head())
    photonField = []
    energy = []
    for i,e in enumerate(eps):
        dens = np.array(list(d.iloc[i])) * 6.2415091e11 / eps[i]**2  # [1/eVcm^3]
        photonField.append(dens)
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Dominguez11(fileDir, outDir):
    name = "IRB_Dominguez11"
    info = "# EBL intensities for the paper >Extragalactic background light inferred from AEGIS galaxy-SED-type fractions<, A. Dominguez et al., 2011, MNRAS, 410, 2556"
    redshift = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.9])
    filePath = fileDir + "EBL_Dominguez_2011/ebl_dominguez11.out"
    d = np.genfromtxt(filePath, unpack=False)
    eps = [1.239842/fieldSlice[0]*crp.eV for fieldSlice in d]  # [eV] | 1.238842 = h*c/1µ eV
    #                          nW->W : J->eV : 1/m^3->1/cm^3 : 1/sm^2sr->1/m^3
    n = np.array([fieldSlice[1:] * 1e-9 * crp.eV / 1e6 / eps[i]**2 * (4 * np.pi / crp.c_light) for i, fieldSlice in enumerate(d)])  # [1/eVcm^3] 
    energy = []
    for i,x in enumerate(n):
        energy.append(eps[i]/crp.eV)
    photonField = [x for x in reversed(n)]
    energy = [e for e in reversed(energy)]
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Stecker16_lower(fileDir, outDir):
    name = "IRB_Stecker16_lower"
    info = "# Extragalactic background light model from Stecker et al. 2016, DOI:10.3847/0004-637X/827/1/6 <An Empirical Determination of the Intergalactic Background Light from UV to FIR Wavelengths Using FIR Deep Galaxy Surveys and the Gamma-ray Opacity of the Universe>"
    redshift = np.linspace(0, 5, 501)
    filePath = fileDir + "EBL_Stecker_2016/comoving_enerdens_lo.csv"
    d = np.genfromtxt(filePath, delimiter=',')
    eps = 10**np.arange(-2.84, 1.14001, 0.01)  # [eV]
    nu = eps * crp.eV / crp.c_light  # [Hz]
    photonField = []
    energy = []
    for i,dens in enumerate(d):
        d[i] = d[i] * 6.2415091e11 * nu[i] / eps[i]**2  # 1/eVcm^3
        photonField.append(d[i])
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Stecker16_upper(fileDir, outDir):
    name = "IRB_Stecker16_upper"
    info = "# Extragalactic background light model from Stecker et al. 2016, DOI:10.3847/0004-637X/827/1/6 <An Empirical Determination of the Intergalactic Background Light from UV to FIR Wavelengths Using FIR Deep Galaxy Surveys and the Gamma-ray Opacity of the Universe>"
    redshift = np.linspace(0, 5, 501)
    filePath = fileDir + "EBL_Stecker_2016/comoving_enerdens_up.csv"
    d = np.genfromtxt(filePath, delimiter=',')
    eps = 10**np.arange(-2.84, 1.14001, 0.01)  # [eV]
    nu = eps * crp.eV / crp.c_light  # [Hz]
    photonField = []
    energy = []
    for i,dens in enumerate(d):
        d[i] = d[i] * 6.2415091e11 * nu[i] / eps[i]**2  # 1/eVcm^3
        photonField.append(d[i])
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Kneiske04(fileDir, outDir):
    name = "IRB_Kneiske04"
    info = "# IRO spectrum from Tanja Kneiske et al. obtained from O. E. Kalashev (http://hecr.inr.ac.ru/)"
    redshift = np.linspace(0,5,51)
    filePath = fileDir + "EBL_Kneiske_2004/all_z"
    d = np.genfromtxt(filePath)
    photonField = []
    energy = []
    for i,fieldSlice in enumerate(d):
        for j,entry in enumerate(fieldSlice):
            if j != 0:
                d[i][j] /= 1e6  # 1/eVm^3 -> 1/eVcm^3
        photonField.append(fieldSlice[1:])
        energy.append(fieldSlice[0])
    createField(name, info, energy, redshift, photonField, outDir)


def IRB_Franceschini08(fileDir, outDir):
    name = "IRB_Franceschini08"
    info = "# Extragalactic background light model from Franceschini et al. 2008, DOI:10.1051/0004-6361:200809691, arxiv/0805.1841v2, tables 1 and 2"
    redshift = np.linspace(0., 2., 11)
    fileDir = fileDir + "EBL_Franceschini_2008/"
    fileList = os.listdir(fileDir)
    fieldData = {}
    # read in data. Note that there is a different energy range for each redshift!
    for i, file in enumerate(sorted(fileList)):
        if "README.txt" not in file:
            data = np.genfromtxt(fileDir + file, unpack=True)
            eps = data[0]  # [log10eV]
            dens = data[1]  # [log10(eps*dn/deps in cm^-3)]
            fieldData[np.round(redshift[i], 2)] = (eps, dens)

    # determine global min and max eps and then generate a logspace from them.
    epsMin, epsMax = [], []
    for key in fieldData.keys():
        epsMin.append(min(fieldData[key][0]))
        epsMax.append(max(fieldData[key][0]))
    energy = 10 ** np.linspace(min(epsMin), max(epsMax), 41)  # [eV]

    # interpolate all redshift-dependent photon field densities on this range
    photonField = []
    for key in fieldData.keys():
        dens = np.interp(energy, 10 ** fieldData[key][0], fieldData[key][1])
        dens = 10 ** dens / energy  # log10(e*dn/de in 1/cm^3) -> 1/(eVcm^3)
        photonField.append(dens)
    photonField = np.swapaxes(np.array(photonField), 0, 1)
    photonField = [zSlice / (1 + redshift)**3 for zSlice in photonField]  # make comoving
    createField(name, info, energy, redshift, photonField, outDir)


def CMB(outDir):
    name = "CMB"
    info = "# Cosmic Microwave Background, T_CMB = 2.72548 K"
    redshift = None
    eps = np.logspace(-10, -1, 101) * crp.eV
    T_CMB = 2.72548
    dnde = lambda e: 8 * np.pi / crp.c_light**3 / crp.h_planck**3 * e**2 / (np.exp(e/(crp.k_boltzmann*T_CMB)) - 1)
    photonField = [[d * crp.eV * cm3] for d in dnde(eps)]  # [1/eVcm^3]
    energy = eps / crp.eV  # [eV]
    createField(name, info, energy, redshift, photonField, outDir)


def URB_Protheroe96(outDir):
    name = "URB_Protheroe96"
    info = "# Universal Radio Background from Protheroe & Bierman Astroparticle Physics 6 (1996) 45."
    redshift = None
    PB = pf.URB_Protheroe96()
    logEmin = np.log10(PB.getEmin())
    logEmax = np.log10(PB.getEmax())
    eps = np.logspace(logEmin, logEmax, 101)
    photonField = np.array([PB.getDensity(e) for e in eps]) / (cm3**-1  * crp.eV**-1)
    energy = eps / crp.eV
    createField(name, info, energy, redshift, photonField, outDir)


def URB_Fixsen11(outDir):
    name = "URB_Fixsen11"
    info = "# Universal Radio Background as measured by ARCADE2 D. J. Fixsen et al. ApJ 734 (2011) 5"
    redshift = None
    FX = pf.URB_Fixsen11()
    logEmin = np.log10(FX.getEmin())
    logEmax = np.log10(FX.getEmax())
    eps = np.logspace(logEmin, logEmax, 101)
    photonField = np.array([FX.getDensity(e) for e in eps]) / (cm3**-1 * crp.eV**-1)
    energy = eps / crp.eV
    createField(name, info, energy, redshift, photonField, outDir)


def URB_Nitu21(outDir):
    name = "URB_Nitu21"
    info = "# Universal Radio Background from Nitu et al. Astropart. Phys. 126 (2021) 102532"
    redshift = None
    NT = pf.URB_Nitu21()
    logEmin = np.log10(NT.getEmin())
    logEmax = np.log10(NT.getEmax())
    eps = np.logspace(logEmin, logEmax, 101)
    photonField = np.array([NT.getDensity(e) for e in eps]) / (cm3**-1  * crp.eV**-1)
    energy = eps / crp.eV
    createField(name, info, energy, redshift, photonField, outDir)


def createField(name, info, energy, redshift, photonDensity, outDir):
    try:
        git_hash = gh.get_git_revision_hash()
        addHash = True
    except:
        addHash = False

    with open(outDir + "/" + name + "_photonEnergy.txt", 'w') as f:
        f.write(info+'\n')
        if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
        f.write("# photon energies in [J]\n")
        for e in energy:
            f.write("{}\n".format(e * crp.eV))  # [J]
    if redshift is not None:
        with open(outDir + "/" + name + "_redshift.txt", 'w') as f:
            f.write(info+'\n')
            if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
            f.write("# redshift\n")
            for z in redshift:
                f.write("{}\n".format(np.round(z, 2)))
    with open(outDir + "/" + name + "_photonDensity.txt", 'w') as f:
        f.write(info+'\n')
        if addHash: f.write("# Produced with crpropa-data version: "+git_hash+"\n")
        f.write("# Comoving photon number density in [m^-3], format: d(e1,z1), ... , d(e1,zm), d(e2,z1), ... , d(e2,zm), ... , d(en,zm)\n")
        for i, densSlice in enumerate(photonDensity):
            for d in densSlice:
                f.write("{}\n".format(d * energy[i] / cm3))  # [# / m^3], comoving
    print("done: " + name)


def process():
    """Processing of the listed photon fields
    
    NOTE: Make sure that this list is synchronized
    with the photon fields used in calc_all.py.
    """
    
    outDir = "data/Scaling"
    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    inFileDir = "tables/"
    IRB_Dominguez11(inFileDir, outDir)
    IRB_Finke10(inFileDir, outDir)
    IRB_Franceschini08(inFileDir, outDir)
    IRB_Gilmore12(inFileDir, outDir)
    IRB_Kneiske04(inFileDir, outDir)
    IRB_Stecker05(inFileDir, outDir)
    IRB_Stecker16_lower(inFileDir, outDir)
    IRB_Stecker16_upper(inFileDir, outDir)
    CMB(outDir)
    URB_Protheroe96(outDir)
    URB_Fixsen11(outDir)
    URB_Nitu21(outDir)

if __name__ == "__main__":
    process()
