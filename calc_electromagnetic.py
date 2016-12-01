from __future__ import division
from numpy import *
import interactionRate
import photonField


eV = 1.60217657E-19  # [J]
me2 = (510.998918E3 * eV)**2  # squared electron mass [J^2/c^4]
sigmaThompson = 6.6524E-29  # Thompson cross section [m^2]
alpha = 1 / 137.035999074  # fine structure constant

def sigmaPP(s):
    """ Pair production cross section (Bethe-Heitler), see Lee 1996. """
    smin = 4 * me2
    if (s < smin):
        return 0
    else:
        b = sqrt(1 - smin / s)
        return sigmaThompson * 3/16 * (1 - b**2) * ((3 - b**4) * log((1 + b) / (1 - b)) - 2 * b * (2 - b**2))

def sigmaICS(s):
    """ Inverse Compton scattering cross sections, see Lee 1996. """
    smin = me2
    if (s < smin):
        return 0
    else:
      b = (s - smin) / (s + smin)
      return sigmaThompson * 3/8 * smin / s / b * (2 / b / (1 + b) * (2 + 2 * b - b**2 - 2 * b**3) - (2 - 3 * b**2 - b**3) / b**2 * log((1 + b) / (1 - b)))

def sigmaTPP(s):
    """ Triplet-pair production cross section, see Lee 1996. """
    beta = 28/9 * log(s / me2) - 218/27
    if beta < 0:
        return 0
    else:
        return sigmaThompson * 3/8 / pi * alpha * beta

def sigmaDPP(s):
    """ Double-pair production cross section, see R.W. Brown eq. (4.5) with k^2 = q^2 = 0. """
    smin = 16 * me2
    if (s < smin):
        return 0
    else:
        return 6.45E-34 * (1 - smin / s)**6

"""
# ----------------------------------------------------------------
# Interaction rates
# ----------------------------------------------------------------
print ('Calculate interaction rates')

def saveRate(E, sigma, skin, field, name):
    if name in ('EMInverseComptonScattering', 'EMTripletPairProduction'):
        xs = array([sigma(s) for s in skin + me2])
    else:
        xs = array([sigma(s) for s in skin])

    rate = interactionRate.calc_rate_s(skin, xs, E, field)

    savetxt('data/%s_%s.txt' % (name, field.name), 
        c_[log10(E / eV), rate], 
        fmt = '%.1f\t%.6g', 
        header = 'log10(E/eV)\t1/lambda [1/Mpc]\n%s\n' % field.info)

# tabulated energies for which the interaction rate is calculated
E = logspace(9, 23, 701) * eV

# tabulated values of skin = s - (mc^2)^2, where s is the Mandelstam s
# Note: the integration method requires log-spaced values
skin1 = logspace(6, 23, 2049) * eV**2  # photons
skin2 = logspace(6, log10(1E23 - me2/eV**2), 2049) * eV**2  # electrons

fields = [
    photonField.CMB(),
    photonField.URB_Protheroe96(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()
    ]

for field in fields:
    print (' -- %s ' % field.name)
    saveRate(E, sigmaPP,  skin1, field, 'EMPairProduction')
    saveRate(E, sigmaDPP, skin1, field, 'EMDoublePairProduction')
    saveRate(E, sigmaICS, skin2, field, 'EMInverseComptonScattering')
    saveRate(E, sigmaTPP, skin2, field, 'EMTripletPairProduction')
"""


# ----------------------------------------------------------------
# Cumulative differential interaction rates
# ----------------------------------------------------------------
print ('Calculate cumulative differential interaction rates')

def saveCDF(E, sigma, skin, field, name):
    if name in ('EMInverseComptonScattering', 'EMTripletPairProduction'):
        xs = array([sigma(s) for s in skin + me2])
    else:
        xs = array([sigma(s) for s in skin])
    
    rate = interactionRate.calc_diffrate_s(skin, xs, E, field)
    # lE = repeat(log10(E/eV), len(skin))
    # ls = repeat(log10(skin/eV**2)[newaxis,:], len(E), axis=0).flatten()
    
    savetxt('data/%s_%s.txt' % (name, field.name), 
        c_[lE, ls, rate], 
        fmt = '%.1f\t%g\t%.6g', 
        header = 'log10(E/eV)\tlog10(s_kin/eV^2)\t(1/lambda)_cumulative [1/Mpc]\n%s\n' % field.info)


# tabulated energies and skin values for which the cumulative differential interaction rate is calculated
E = logspace(9, 23, 141) * eV
skin1 = logspace(6, 23, 171) * eV**2 # photons
skin2 = logspace(6, log10(1E23 - me2/eV**2), 171) * eV**2 # electrons

fields = [
    photonField.CMB(),
    # photonField.EBL_Kneiske04(),
    # photonField.EBL_Stecker05(),
    # photonField.EBL_Franceschini08(),
    # photonField.EBL_Finke10(),
    # photonField.EBL_Dominguez11(),
    # photonField.EBL_Gilmore12()
    ]

for field in fields:
    print (' -- %s ' % field.name)
    saveCDF(E, sigmaPP,  skin1, field, 'EMPairProduction_CDF')
    # saveCDF(E, sigmaICS, skin2, field, 'EMInverseComptonScattering_CDF')
    # saveCDF(E, sigmaTPP, skin2, field, 'EMTripletPairProduction_CDF')


# # Consider different s-range for URB
# field = photonField.URB_Protheroe96()
# print (' -- %s ' % field.name)
# skin1 = logspace(6, 17, 111) + log10(eV**2)  # photons
# skin2 = logspace(6, log10(1E17 - me2/eV**2), 111) + log10(eV**2)  # electrons
# saveCDF(E, sigmaPP,  skin1, field, 'EMPairProduction_CDF')
# saveCDF(E, sigmaICS, skin2, field, 'EMInverseComptonScattering_CDF')
# saveCDF(E, sigmaTPP, skin2, field, 'EMTripletPairProduction_CDF')
