from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19
gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


# cross sections for A < 12 (various sources)
ddir1 = 'tables/PD_external/'
isotopes1 = genfromtxt(ddir1 + 'isotopes.txt')
x = genfromtxt(ddir1+'eps.txt') * eV * 1e6  # [J]
n = len(x)
d1sum = genfromtxt(ddir1+'xs_sum.txt',
    dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
d1exc = genfromtxt(ddir1+'xs_excl.txt',
    dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%n)])
eps1 = iR.romb_pad_logspaced(x, 513)  # padding
xs1sum = array([iR.romb_pad_zero(x, 513) for x in d1sum['xs']])*1e-31
xs1exc = array([iR.romb_pad_zero(x, 513) for x in d1exc['xs']])*1e-31


# cross sections for A >= 12 (TALYS)
ddir2 = 'tables/PD_Talys1.6_Khan/'
# ddir2 = 'tables/PD_Talys1.6_Geant4/'
isotopes2 = genfromtxt(ddir2 + 'isotopes.txt')
x = genfromtxt(ddir2+'eps.txt') * eV * 1e6  # [J]
n = len(x)
d2sum = genfromtxt(ddir2+'xs_sum.txt',
    dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
d2exc = genfromtxt(ddir2+'xs_thin.txt',
    dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%n)])
eps2 = iR.romb_pad_logspaced(x, 513)  # padding
xs2sum = array([iR.romb_pad_zero(x, 513) for x in d2sum['xs']])*1e-31
xs2exc = array([iR.romb_pad_zero(x, 513) for x in d2exc['xs']])*1e-31


fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12()]


for field in fields:
    print field.name

    # Calculate total interaction rate
    R1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1sum])
    R2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2sum])

    # save
    fname = 'data/pd_%s.txt' % field.name
    output = r_[
        c_[d1sum['Z'], d1sum['N'], R1],
        c_[d2sum['Z'], d2sum['N'], R2]]
    fmt = '%i\t%i' + '\t%g'*201
    hdr = 'Photo-disintegration with the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)


    # Calculate branching ratios
    # for A < 12
    B1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1exc])
    for (Z, N, A) in isotopes1:
        s = (d1exc['Z'] == Z) * (d1exc['N'] == N)
        B1[s] /= sum(B1[s], axis=0)
    B1[isnan(B1)] = 0  # set to 0 when total cross section is 0

    # for A > 12
    B2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2exc])
    for (Z, N, A) in isotopes2:
        s = (d2exc['Z'] == Z) * (d2exc['N'] == N)
        B2[s] /= sum(B2[s], axis=0)
    B2[isnan(B2)] = 0  # set to 0 when total cross section is 0

    # save
    fname = 'data/pd_branching_%s.txt'%field.name
    output = r_[
        c_[d1exc['Z'], d1exc['N'], d1exc['ch'], B1],
        c_[d2exc['Z'], d2exc['N'], d2exc['ch'], B2]]
    fmt = '%i\t%i\t%06d' + '\t%g'*201
    hdr = 'Photo-disintegration with the %s\nZ, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)
