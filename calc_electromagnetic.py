from numpy import *
import interactionRate
import photonField


eV = 1.60217657e-19  # [J]
me2 = 0.510998918e6**2  # squared electron mass [eV^2/c^4]
sigmaThompson = 6.6524e-29  # Thompson cross section [m^2]
alpha = 1. / 137.035999074  # fine structure constant

def sigmaPP(s):
    """
    Pair production cross section (Bethe-Heitler), see Lee 1996.
    Returns cross sections [1/m^2] for chosen s range [J^2].
    """
    if (s < 4.*me2 * eV**2):
        return 0.
    else:
        beta = sqrt(1 - 4*me2 * eV**2 /s)
        return sigmaThompson * 3./16.*(1-beta**2)*((3-beta**4)*log((1+beta)/(1-beta))-2*beta*(2-beta**2))

def sigmaICS(s):
    """
    Inverse Compton scattering cross sections, see Lee 1996.
    Returns cross sections [1/m^2] for chosen s range [J^2].
    """
    if (s < 1.0 * me2 *eV**2):
        return 0.
    else:
      beta = (s - me2 *eV**2)/(s + me2 *eV**2)
      return sigmaThompson * 3./8. * me2 * eV**2 / s / beta * (2./beta/(1+beta)*(2.+2.*beta-beta**2-2.*beta**3)-1./beta**2 * (2.-3.*beta**2- beta**3)*log((1.+beta)/(1.-beta)))

def sigmaTPP(s):
    """
    Triplet-pair production cross section, see Lee 1996.
    Returns cross sections [1/m^2] for chosen s range [J^2].
    """
    if (28/9*log(s/me2 /eV**2)- 218./27. < 0.):
        return 0.
    else:
        return sigmaThompson *3.*alpha/8./pi*(28/9*log(s/me2 /eV**2)- 218./27.)

def sigmaDPP(s):
    """
    Double-pair production cross section, see R.W. Brown eq. (4.5) with k^2 = q^2 = 0.
    Returns cross sections [1/m^2] for chosen s range [J^2].
    """
    if (s < 16 * me2 * eV **2):
        return 0.
    else:
        # exponent = 1 instead of 6 results in better reproduction of EleCa data
        return 6.45*1e-34 *(1.- 16.*me2 *eV**2 / s)**6


# ----------------------------------------------------------------
# Interaction rates
# ----------------------------------------------------------------
print ('Calculate interaction rates')

def saveRate(E, sigma, skin, field, name):
    s = skin.copy()
    if (name == 'EMInverseComptonScattering' or name == 'EMTripletPairProduction'):
      s += me2 * eV**2
    xs   = array([sigma(si) for si in s])
    rate = interactionRate.calc_rate_s(skin, xs, E, field)
    data = c_[log10(E / eV), rate]
    fname  = 'data/%s_%s.txt' % (name, field.name)
    header = 'log10(E/eV)\t1/lambda [1/Mpc]\n%s\n' % field.info
    fmt = '%.2f\t%.6e'
    savetxt(fname, data, fmt=fmt, header=header)

# Mandelstam s - (mc^2)^2
# important: tabulate skin in logspace and add resp. mass term to obtain s
# because the integration carried out requires logspaced skin values.
# Otherwise (logspace s from mass term to 10**23 eV and subtract mass term
# to obtain skin) skin is not logspaced for small values.
skin1 = logspace(log10(1E7*eV**2), log10( 10E23*eV**2), 1025)  # photons
skin2 = logspace(log10(1E7*eV**2), log10((10E23-me2)*eV**2), 1025)  # electrons
E = logspace(15, 23, 801) * eV  # energy range of interacting particle

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
    print (field.name)
    saveRate(E, sigmaPP,  skin1, field, 'EMPairProduction')
    saveRate(E, sigmaDPP, skin1, field, 'EMDoublePairProduction')
    saveRate(E, sigmaICS, skin2, field, 'EMInverseComptonScattering')
    saveRate(E, sigmaTPP, skin2, field, 'EMTripletPairProduction')


# ----------------------------------------------------------------
# Cumulative differential interaction rates
# ----------------------------------------------------------------
print ('Calculate cumulative differential interaction rates')

def saveCDF(E, sigma, skin, field, name):
    s = skin.copy()
    if name in ('EMInverseComptonScattering_CDF', 'EMTripletPairProduction_CDF'):
        s += me2 * eV**2
    xs   = array([sigma(si) for si in s])
    rate = interactionRate.calc_diffrate_s(skin, xs, E, field)
    lE = repeat(log10(E/eV), len(skin))
    ls = repeat(log10(skin/eV**2)[newaxis,:], len(E), axis=0).flatten()
    data = c_[lE, ls, rate]
    fname  = 'data/%s_%s.txt' % (name, field.name)
    header = 'log10(E/eV)\tlog10(s_kin/eV^2)\t(1/lambda)_cumulative [1/Mpc]\n%s\n' % field.info
    fmt = '%.2f\t%.5e\t%.5e'
    savetxt(fname, data, fmt=fmt, header=header)

skin1 = logspace(log10(1E7*eV**2), log10( 1E23*eV**2), 500)
skin2 = logspace(log10(1E7*eV**2), log10((1E23-me2)*eV**2), 500)
E = logspace(15, 23, 81) * eV  # energy range of interacting particle

fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()
    ]

for field in fields:
    print (field.name)
    saveCDF(E, sigmaPP,  skin1, field, 'EMPairProduction_CDF')
    saveCDF(E, sigmaICS, skin2, field, 'EMInverseComptonScattering_CDF')
    saveCDF(E, sigmaTPP, skin2, field, 'EMTripletPairProduction_CDF')


# Consider different s-range for URB
field = photonField.URB_Protheroe96()
skin1 = logspace(log10(1E7*eV**2), log10( 1E17 * eV**2), 500)
skin2 = logspace(log10(1E7*eV**2), log10((1E17 - me2)*eV**2), 500)
saveCDF(E, sigmaPP,  skin1, field, 'EMPairProduction_CDF')
saveCDF(E, sigmaICS, skin2, field, 'EMInverseComptonScattering_CDF')
saveCDF(E, sigmaTPP, skin2, field, 'EMTripletPairProduction_CDF')
