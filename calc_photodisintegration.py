from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19
gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


# Load cross sections for A < 12
ddir1 = 'tables/PD_external/'
isotopes1 = genfromtxt(ddir1 + 'isotopes.txt')
eps = genfromtxt(ddir1 + 'eps.txt')
d1sum = genfromtxt(ddir1 + 'xs_sum.txt', dtype=[('Z',int), ('N',int), ('xs','%if8'%len(eps))])
d1exc = genfromtxt(ddir1 + 'xs_excl.txt', dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%len(eps))])
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps1 = iR.romb_pad_logspaced(eps, 513) * eV * 1e6
xs1sum = array([iR.romb_pad_zero(x, 513) for x in d1sum['xs']]) * 1e-31
xs1exc = array([iR.romb_pad_zero(x, 513) for x in d1exc['xs']]) * 1e-31


# Load cross sections for A >= 12 (TALYS)
ddir2 = 'tables/PD_Talys1.8_Khan/'
isotopes2 = genfromtxt(ddir2 + 'isotopes.txt')
eps = genfromtxt(ddir2+'eps.txt')
d2sum = genfromtxt(ddir2+'xs_pd_sum.txt', dtype=[('Z',int), ('N',int), ('xs','%if8'%len(eps))])
d2exc = genfromtxt(ddir2+'xs_pd_thin.txt', dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%len(eps))])
# Only consider cross sections for A > 12
d2sum = d2sum[(d2sum['Z'] + d2sum['N']) >= 12]
d2exc = d2exc[(d2exc['Z'] + d2exc['N']) >= 12]
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps2 = iR.romb_pad_logspaced(eps, 513) * eV * 1e6
xs2sum = array([iR.romb_pad_zero(x, 513) for x in d2sum['xs']]) * 1e-31
xs2exc = array([iR.romb_pad_zero(x, 513) for x in d2exc['xs']]) * 1e-31


# Load cross sections with photon emission
d3sum = genfromtxt(ddir2+'xs_photon_sum.txt', dtype=[('Z',int), ('N',int), ('Zd',int), ('Nd',int), ('xs','%if8'%len(eps))])
d3exc = genfromtxt(ddir2+'xs_photon_thin.txt', dtype=[('Z',int), ('N',int), ('Zd',int), ('Nd',int), ('Ephoton',float), ('xs','%if8'%len(eps))])
# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration and convert to SI units
eps3 = eps2
xs3sum = array([iR.romb_pad_zero(x, 513) for x in d3sum['xs']]) * 1e-31
xs3exc = array([iR.romb_pad_zero(x, 513) for x in d3exc['xs']]) * 1e-31



fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()]


for field in fields:
    print (field.name)

    # Calculate total interaction rates
    R1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1sum])
    R2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2sum])

    savetxt('data/pd_%s.txt' % field.name,
        r_[c_[d1sum['Z'], d1sum['N'], R1], c_[d2sum['Z'], d2sum['N'], R2]],
        fmt='%i\t%i' + '\t%g'*201,
        header='Photodisintegration by the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info)


    # Calculate branching ratios from exclusive interaction rates
    B1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1exc])
    B2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2exc])
    for (Z, N, A) in isotopes1:
        s = (d1exc['Z'] == Z) * (d1exc['N'] == N)
        B1[s] /= sum(B1[s], axis=0)
    for (Z, N, A) in isotopes2:
        s = (d2exc['Z'] == Z) * (d2exc['N'] == N)
        B2[s] /= sum(B2[s], axis=0)
    B1 = nan_to_num(B1)  # set to 0 when total cross section is 0
    B2 = nan_to_num(B2)

    savetxt('data/pd_branching_%s.txt'%field.name,
        r_[c_[d1exc['Z'], d1exc['N'], d1exc['ch'], B1], c_[d2exc['Z'], d2exc['N'], d2exc['ch'], B2]],
        fmt='%i\t%i\t%06d' + '\t%g'*201,
        header='Photo-disintegration with the %s\nZ, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps' % field.info)


    # Calculate photon emission probabilities
    R3 = array([iR.invMFP_fast(eps3, x, gamma, field) for x in xs3sum])
    B3 = array([iR.invMFP_fast(eps3, x, gamma, field) for x in xs3exc])
    for i in range(len(d3sum)):
        s = (d3exc['Z'] == d3sum['Z'][i]) * (d3exc['N'] == d3sum['N'][i]) * (d3exc['Zd'] == d3sum['Zd'][i]) * (d3exc['Nd'] == d3sum['Nd'][i])
        B3[s] /= R3[i]
    B3 = nan_to_num(B3)

    savetxt('data/pd_photon_emission_%s.txt' % field.name,
        c_[d3exc['Z'], d3exc['N'], d3exc['Zd'], d3exc['Nd'], d3exc['Ephoton']*1e6, B3],
        fmt='%i\t%i\t%i\t%i\t%g' + '\t%g'*201,
        header='Emission probabilities of photons with discrete energies via photo-disintegration with the %s\nZ, N, Z_daughter, N_daughter, Ephoton [eV], emission probability for log10(gamma) = 6-14 in 201 steps' % field.info)
