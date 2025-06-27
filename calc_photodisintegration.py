import numpy as np
import os
import interactionRate
import gitHelp as gh
from calc_all import fields_cmbebl, fields_urb, reduced_fields
from units import eV

cdir = os.path.split(__file__)[0]

gamma = np.logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


# ----------------------------------------------------
# Load cross sections for A < 12
# ----------------------------------------------------
ddir1 = os.path.join(cdir, 'tables/PD_external/')
isotopes1 = np.genfromtxt(ddir1 + 'isotopes.txt')
eps = np.genfromtxt(ddir1 + 'eps.txt')
d1sum = np.genfromtxt(ddir1 + 'xs_sum.txt', dtype=[('Z', int), ('N', int), ('xs', '%if8' % len(eps))])
d1exc = np.genfromtxt(ddir1 + 'xs_excl.txt', dtype=[('Z', int), ('N', int), ('ch', int), ('xs', '%if8' % len(eps))])
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps1 = interactionRate.romb_pad_logspaced(eps, 513) * eV * 1e6
xs1sum = np.array([interactionRate.romb_pad_zero(x, 513) for x in d1sum['xs']]) * 1e-31
xs1exc = np.array([interactionRate.romb_pad_zero(x, 513) for x in d1exc['xs']]) * 1e-31


# ----------------------------------------------------
# Load cross sections for A >= 12 (TALYS)
# ----------------------------------------------------
ddir2 = os.path.join(cdir, 'tables/PD_Talys1.8_Khan/')
isotopes2 = np.genfromtxt(ddir2 + 'isotopes.txt')
eps = np.genfromtxt(ddir2 + 'eps.txt')
d2sum = np.genfromtxt(ddir2 + 'xs_pd_sum.txt', dtype=[('Z', int), ('N', int), ('xs', '%if8' % len(eps))])
d2exc = np.genfromtxt(ddir2 + 'xs_pd_thin.txt', dtype=[('Z', int), ('N', int), ('ch', int), ('xs', '%if8' % len(eps))])
# Only consider cross sections for A > 12
d2sum = d2sum[(d2sum['Z'] + d2sum['N']) >= 12]
d2exc = d2exc[(d2exc['Z'] + d2exc['N']) >= 12]
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps2 = interactionRate.romb_pad_logspaced(eps, 513) * eV * 1e6
xs2sum = np.array([interactionRate.romb_pad_zero(x, 513) for x in d2sum['xs']]) * 1e-31
xs2exc = np.array([interactionRate.romb_pad_zero(x, 513) for x in d2exc['xs']]) * 1e-31


# ----------------------------------------------------
# Load cross sections for A >= 56 (TALYS)
# ----------------------------------------------------
ddir4 = os.path.join(cdir, 'tables/PD_Talys1.9/')
isotopes4 = np.genfromtxt(ddir4 + 'isotopes.txt')
eps4 = np.genfromtxt(ddir4 + 'eps.txt')
d4sum = np.genfromtxt(ddir4 + 'xs_pd_sum.txt', dtype=[('Z', int), ('N', int), ('xs', '%if8' % len(eps))])
d4exc = np.genfromtxt(ddir4 + 'xs_pd_thin.txt', dtype=[('Z', int), ('N', int), ('ch', int), ('xs', '%if8' % len(eps))])
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps4 = interactionRate.romb_pad_logspaced(eps4, 513) * eV * 1e6
xs4sum = np.array([interactionRate.romb_pad_zero(x, 513) for x in d4sum['xs']]) * 1e-31
xs4exc = np.array([interactionRate.romb_pad_zero(x, 513) for x in d4exc['xs']]) * 1e-31


# ----------------------------------------------------
# Load cross sections with photon emission
# ----------------------------------------------------
d3sum = np.genfromtxt(ddir2 + 'xs_photon_sum.txt', dtype=[('Z', int), ('N', int), ('Zd', int), ('Nd', int), ('xs', '%if8' % len(eps))])
d3exc = np.genfromtxt(ddir2 + 'xs_photon_thin.txt', dtype=[('Z', int), ('N', int), ('Zd', int), ('Nd', int), ('Ephoton', float), ('xs', '%if8' % len(eps))])
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps3 = eps2
xs3sum = np.array([interactionRate.romb_pad_zero(x, 513) for x in d3sum['xs']]) * 1e-31
xs3exc = np.array([interactionRate.romb_pad_zero(x, 513) for x in d3exc['xs']]) * 1e-31


# ----------------------------------------------------
# Calculate interaction rates and branching ratios
# ----------------------------------------------------
def processRate(field):

    # output folder
    folder = 'data/Photodisintegration'
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Calculate total interaction rates
    R1 = np.array([interactionRate.calc_rate_eps(eps1, x, gamma, field) for x in xs1sum])
    R2 = np.array([interactionRate.calc_rate_eps(eps2, x, gamma, field) for x in xs2sum])
    R4 = np.array([interactionRate.calc_rate_eps(eps4, x, gamma, field) for x in xs4sum])
    
    try:
        git_hash = gh.get_git_revision_hash()
        header = ("Photodisintegration rate with the %s\n"% field.info
                  +"Produced with crpropa-data version: "+git_hash+"\n"
                  +"Z, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps" )
    except:
        header = ("Photodisintegration rate with the %s\n"
                  "Z, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps" % field.info)
    np.savetxt(
        folder + '/rate_%s.txt' % field.name,
        np.r_[np.c_[d1sum['Z'], d1sum['N'], R1], np.c_[d2sum['Z'], d2sum['N'], R2]],
        fmt='%i\t%i' + '\t%.9g' * 201,
        header=header)
    
    np.savetxt(
        folder + '/rate_%s_superheavy.txt' % field.name,
        np.r_[np.c_[d4sum['Z'], d4sum['N'], R4]],
        fmt='%i\t%i' + '\t%.9g' * 201,
        header=header)

    # Calculate branching ratios from exclusive interaction rates
    B1 = np.array([interactionRate.calc_rate_eps(eps1, x, gamma, field) for x in xs1exc])
    B2 = np.array([interactionRate.calc_rate_eps(eps2, x, gamma, field) for x in xs2exc])
    B4 = np.array([interactionRate.calc_rate_eps(eps4, x, gamma, field) for x in xs4exc])
    for (Z, N, A) in isotopes1:
        s = (d1exc['Z'] == Z) * (d1exc['N'] == N)
        B1[s] /= np.sum(B1[s], axis=0)
    for (Z, N, A) in isotopes2:
        s = (d2exc['Z'] == Z) * (d2exc['N'] == N)
        B2[s] /= np.sum(B2[s], axis=0)
    for (Z, N, A) in isotopes4:
        s = (d4exc['Z'] == Z) * (d4exc['N'] == N)
        B4[s] /= np.sum(B4[s], axis=0)
    B1 = np.nan_to_num(B1)  # set to 0 when total cross section is 0
    B2 = np.nan_to_num(B2)
    B4 = np.nan_to_num(B4)
    
    try:
        git_hash = gh.get_git_revision_hash()
        header = ("Photodisintegration rate with the %s\n"% field.info
                  +"Produced with crpropa-data version: "+git_hash+"\n"
                  +"Z, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps" )
    except:
        header = ("Photodisintegration rate with the %s\n"
                  "Z, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps" % field.info)
    np.savetxt(
        folder + '/branching_%s.txt' % field.name,
        np.r_[np.c_[d1exc['Z'], d1exc['N'], d1exc['ch'], B1], np.c_[d2exc['Z'], d2exc['N'], d2exc['ch'], B2]],
        fmt='%i\t%i\t%06d' + '\t%.9g' * 201,
        header=header)

    np.savetxt(
        folder + '/branching_%s_superheavy.txt' % field.name,
        np.r_[np.c_[d4exc['Z'], d4exc['N'], d4exc['ch'], B4]],
        fmt='%i\t%i\t%06d' + '\t%.9g' * 201,
        header=header)


# ----------------------------------------------------
# Calculate photon emission probabilities
# ----------------------------------------------------
def processEmission(field):
    """ calculate photon emission probabilities """

    folder = "data/Photodisintegration/"
    if not os.path.isdir:
        os.makedirs(folder)

    R3 = np.array([interactionRate.calc_rate_eps(eps3, x, gamma, field) for x in xs3sum])
    B3 = np.array([interactionRate.calc_rate_eps(eps3, x, gamma, field) for x in xs3exc])
    for i in range(len(d3sum)):
        s = (d3exc['Z'] == d3sum['Z'][i]) * (d3exc['N'] == d3sum['N'][i]) * (d3exc['Zd'] == d3sum['Zd'][i]) * (d3exc['Nd'] == d3sum['Nd'][i])
        B3[s] /= R3[i]
    B3 = np.nan_to_num(B3)
    
    try:
        git_hash = gh.get_git_revision_hash()
        header = ("Emission probabilities of photons with discrete energies via photo-disintegration with the%s\n"% field.info
                  +"Produced with crpropa-data version: "+git_hash+"\n"
                  +"Z, N, Z_daughter, N_daughter, Ephoton [eV], emission probability for log10(gamma) = 6-14 in 201 steps" )
    except:
        header = ("Emission probabilities of photons with discrete energies via photo-disintegration with the%s\n"
                  "Z, N, Z_daughter, N_daughter, Ephoton [eV], emission probability for log10(gamma) = 6-14 in 201 steps" % field.info)


    np.savetxt(
        folder + 'photon_emission_%s.txt' % field.name[:3],
        np.c_[d3exc['Z'], d3exc['N'], d3exc['Zd'], d3exc['Nd'], d3exc['Ephoton'] * 1e6, B3],
        fmt='%i\t%i\t%i\t%i\t%.9g' + '\t%.9g' * 201,
        header=header)

if __name__ == "__main__":

    for field in fields_cmbebl+fields_urb:
        print('{} interaction rate'.format(field.name))
        processRate(field)

    # only representive fields of each type (CMB, EBL, URB) to save data space
    for field in  reduced_fields:
        print('{} secondary yields'.format(field.name))
        processEmission(field)

