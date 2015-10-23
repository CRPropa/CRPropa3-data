from numpy import *
import interactionRate as iR
import photonField


eV = 1.60217657e-19

# proton / neutron cross sections [1/m^2] for tabulated energies [J]
# truncate to largest length 2^i + 1 for Romberg integration

# photo-pair
eps1 =   # [J]
xs1  =   # [m^2]

# double photo pair
eps2 =
xs2  =

# inverse compton
eps3 =
xs3  =

# triplet pair
eps4 =
xs4  =


lgamma = linspace(6, 16, 251)
gamma  = 10**lgamma

fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()]

for field in fields:
    r1 = iR.invMFP_fast(eps1[0:2049], xs1[0:2049], gamma, field)
    r2 = iR.invMFP_fast(eps2[0:2049], xs2[0:2049], gamma, field)

    fname = 'data/ppp_%s.txt' % field.name
    data  = c_[lgamma, r1, r2]
    fmt   = '%.2f\t%.6e\t%.6e'
    header = 'Photo-pion interaction rate with the %s\nlog10(gamma)\t1/lambda_proton [1/Mpc]\t1/lambda_neutron [1/Mpc]'%field.info
    savetxt(fname, data, fmt=fmt, header=header)
