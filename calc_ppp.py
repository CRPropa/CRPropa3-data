from numpy import *
import interactionRate as iR
import photonField


eV = 1.60217657e-19

# proton / neutron cross sections [1/m^2] for tabulated energies [J]
# truncate to largest length 2^i + 1 for Romberg integration
d = genfromtxt('tables/PPP/xs_proton.txt', unpack=True)
eps1 = d[0] * 1e9 * eV  # [J]
xs1  = d[1] * 1e-34  # [m^2]

d = genfromtxt('tables/PPP/xs_neutron.txt', unpack=True)
eps2 = d[0] * 1e9 * eV  # [J]
xs2  = d[1] * 1e-34  # [m^2]

lgamma = linspace(6, 16, 251)
gamma  = 10**lgamma

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
    r1 = iR.invMFP_fast(eps1[0:2049], xs1[0:2049], gamma, field)
    r2 = iR.invMFP_fast(eps2[0:2049], xs2[0:2049], gamma, field)

    fname = 'data/ppp_%s.txt' % field.name
    data  = c_[lgamma, r1, r2]
    fmt   = '%.2f\t%.6e\t%.6e'
    header = ("Photo-pion interaction rate with the %s\nlog10(gamma)"
              "\t1/lambda_proton [1/Mpc]\t1/lambda_neutron [1/Mpc]"%field.info)
    savetxt(fname, data, fmt=fmt, header=header)

# redshift dependent
for field in fields:

    if field.redshift is None:
        continue

    data = []
    for z in field.redshift:
        r1 = iR.invMFP_fast(eps1[0:2049], xs1[0:2049], gamma, field)
        r2 = iR.invMFP_fast(eps2[0:2049], xs2[0:2049], gamma, field)
        data.append( c_[[z]*len(lgamma), lgamma, r1, r2] )

    data = concatenate( [d for d in data], axis=0 )
    nan_to_num( data )
    fname = 'data/ppp_%sz.txt' % field.name
    fmt   = '%.2f\t%.2f\t%.6e\t%.6e'
    header = ("Photo-pion interaction rate for different redshifts %s\n"
              "z\tlog10(gamma)\t1/lambda_proton [1/Mpc]\t1/lambda_neutron [1/Mpc]"%field.info)
    savetxt(fname, data, fmt=fmt, header=header)
