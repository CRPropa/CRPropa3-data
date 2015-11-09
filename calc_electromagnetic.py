from numpy import *
import interactionRate as iR
import photonField
import matplotlib.pyplot as plt


eV = 1.60217657e-19
ElectronMass = 0.510998918e6 # [eV/c^2]
K_boltz = 8.617342294984e-5  # [eV/K ] Boltzman constant
C_speed = 299792458          # [m/s] speed of light
SigmaThompson = 6.6524e-29   # m**2
alpha = 1. / 137.035999074

# cross sections [1/m^2] for choosen s range [J^2]

def Sigma_PP(s):
    """
    comp. Lee 96
    """
    if (s < 4.*ElectronMass**2 *eV**2):
        return 0.
    else:
        beta = sqrt(1.-4.*ElectronMass**2 *eV**2 /s)
        return SigmaThompson * 3./16.*(1-beta**2)*((3-beta**4)*log((1+beta)/(1-beta))-2*beta*(2-beta**2))

def Sigma_ICS(s):
    """
    comp. Lee 96
    """
    if (s < 1.01 * ElectronMass**2 *eV**2):
        return 0.
    else:
      beta = (s - ElectronMass**2 *eV**2)/(s + ElectronMass**2 *eV**2)
      return SigmaThompson * 3./8. * ElectronMass**2 * eV**2 / s / beta * (2./beta/(1+beta)*(2.+2.*beta-beta**2-2.*beta**3)-1./beta**2 * (2.-3.*beta**2- beta**3)*log((1.+beta)/(1.-beta)))

def Sigma_TPP(s):
    """
    comp. Lee 96
    """
    if (28/9*log(s/ElectronMass**2 /eV**2)- 218./27. < 0.):
        return 0.
    else:
        return SigmaThompson *3.*alpha/8./pi*(28/9*log(s/ElectronMass**2 /eV**2)- 218./27.)

def Sigma_DPP(s):
    """
    comp. R.W. Brown eq. (4.5) k**2 = q**2 = 0 case    
    """
    if (s < 16 * ElectronMass**2 * eV **2):
        return 0.
    else: 
        return 6.45*1e-34 *(1.- 16.*ElectronMass**2 *eV**2 / s)**6
 
# truncate to largest length 2^i + 1 for Romberg integration
s = 10 ** linspace(log10(1.1*ElectronMass**2*eV**2),log10(10**23 * eV**2),1025)
skin_PP_DPP = s
skin_ICS_TPP = s - ElectronMass**2 * eV**2

# choose energy range of interacting particle for which the interaction rate should be tabulated
lEnergy = linspace(12.75, 23, 1026)
Energy  = 10**lEnergy * eV

# photo pair production PP
xs1 = zeros(len(s))
for i in range(0,len(s)):
    xs1[i] = Sigma_PP(s[i])

# inverse compton scattering ICS
xs2 = zeros(len(s))
for i in range(0,len(s)):
  xs2[i] = Sigma_ICS(s[i])  #resulting interaction lengths differ from eleca.dat -> further cross check needed  

# photo double pair production DPP
xs3 = zeros(len(s))
for i in range(0,len(s)):
  xs3[i] = Sigma_DPP(s[i])

# triplet pair production TPP
xs4 = zeros(len(s))
for i in range(0,len(s)):
  xs4[i] = Sigma_TPP(s[i])

fields = [
    photonField.URB_Protheroe96(),
    photonField.CMB(),
    photonField.CRB_Biermann96(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()
    ]

for field in fields:

    r1 = iR.rate(skin_PP_DPP, xs1, Energy, field)
    r2 = iR.rate(skin_ICS_TPP, xs2, Energy, field)
    r3 = iR.rate(skin_PP_DPP, xs3, Energy, field)
    r4 = iR.rate(skin_ICS_TPP, xs4, Energy, field)

    fname = 'data/EleCa/interactionlength_%s.txt' % field.name
    data  = c_[log10(Energy/eV),r1,r2,r3,r4]
    fmt   = '%.3f\t%.6e\t%.6e\t%.6e\t%.6e'
    header = 'Interaction rate with the %s\nlog10(E/eV)\t1/lambda_PP [1/Mpc]\t1/lambda_ICS [1/Mpc]\t1/lambda_DPP [1/Mpc]\t1/lambda_TPP [1/Mpc]'%field.info
    savetxt(fname, data, fmt=fmt, header=header)
