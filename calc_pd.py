from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19
gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


# cross sections for A < 12 (TALYS)
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
idx = d2sum['Z'] + d2sum['N'] >= 12
d2sum = d2sum[idx]
d2exc = genfromtxt(ddir2+'xs_thin.txt',
    dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%n)])
idx = d2exc['Z'] + d2exc['N'] >= 12
d2exc = d2exc[idx]
eps2 = iR.romb_pad_logspaced(x, 513)  # padding
xs2sum = array([iR.romb_pad_zero(x, 513) for x in d2sum['xs']])*1e-31
xs2exc = array([iR.romb_pad_zero(x, 513) for x in d2exc['xs']])*1e-31


# cross section of photon emission after PD for A <= 56 (TALYS)
d_exc_photon = genfromtxt(ddir2+'xs_photon_norm.txt',
    dtype=[('Z',int), ('N',int), ('Zd',int), ('Nd',int), ('Ephoton',double), ('xs','%if8'%n)])
d_sum_photon = genfromtxt(ddir2+'xs_sum_daughter.txt',
    dtype=[('Z',int), ('N',int), ('Zd',int), ('Nd',int), ('xs','%if8'%n)])
xs_sum_photon = array([iR.romb_pad_zero(x, 513) for x in d_sum_photon['xs']])*1e-31
xs_exc_photon = array([iR.romb_pad_zero(x, 513) for x in d_exc_photon['xs']])*1e-31


# cross section of elastic scattering for A <= 56 (TALYS)
x = genfromtxt(ddir2+'eps_elastic.txt') * eV * 1e6  # [J]
n = len(x)
eps_elastic = iR.romb_pad_logspaced(x, 513)  # padding
d_elastic = genfromtxt(ddir2+'xs_elastic.txt',
    dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
xs_elastic = array([iR.romb_pad_zero(x, 513) for x in d_elastic['xs']])*1e-31

# set some isotopes list for use in calculation of interaction rates
iso = vstack(tuple(r) for r in transpose((d_elastic['Z'],d_elastic['N'],d_elastic['Z']+d_elastic['N'])))
idA = iso[:,0] + iso[:,1] < 12
Niso = len(iso)
NisoA = len(iso[idA])
all_isotopes = vstack((isotopes1,isotopes2))
idMiss = zeros(len(all_isotopes),dtype=bool)
for r in iso:
	h = all_isotopes == r
	l = h[:,0]*h[:,1]
	idMiss += l
idMiss = ~idMiss
missing = all_isotopes[idMiss] 
all_isotopes = all_isotopes[~idMiss]


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
    R3 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs_sum_photon])
    R4 = array([iR.invMFP_fast(eps_elastic, x, gamma, field) for x in xs_elastic])

    # Calculate total interaction rate for elastic scattering via TRK scaling of rates for A >= 12
    M = outer(zeros(Niso - NisoA),zeros(len(gamma)))
    m = zeros(len(gamma))

    for i in range(NisoA, Niso):
        M[i - NisoA] = R4[i] / d_elastic['Z'][i] / d_elastic['N'][i] * (d_elastic['Z'][i] + d_elastic['N'][i])
    for i in range(0,len(gamma)):
        m[i] = mean(M[NisoA:,i])
    # use TRK scaling to get interaction length for A < 12 and A >= 12
    for i, (Z,N,A) in enumerate(all_isotopes):
        R4[i] = m * Z * N / A
    # use TRK scaling to get interaction length for A <= 5 (can not be calculated with TALYS)
    R5 = zeros((len(missing),len(R4[0,:])))
    for i, (Z,N,A) in enumerate(missing):
        R5[i] = m * Z * N / A
    # add missing isotopes
    R4 = vstack((R5,R4))
    d_elasticZ, d_elasticN = transpose(vstack((missing,transpose(vstack((d_elastic['Z'],d_elastic['N'],d_elastic['Z']+d_elastic['N'])))))[:,:2])

    # save interaction rate for photo disintegration
    fname = 'data/pd_%s.txt' % field.name
    output = r_[
        c_[d1sum['Z'], d1sum['N'], R1],
        c_[d2sum['Z'], d2sum['N'], R2]]
    fmt = '%i\t%i' + '\t%g'*201
    hdr = 'Photo-disintegration with the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)

    # save interaction rate for elastic scattering with TRK scaling
    fname = 'data/ElasticScattering_%s.txt' % field.name
    output = c_[d_elasticZ, d_elasticN, R4]
    fmt = '%i\t%i' + '\t%g'*201
    hdr = 'Elastic scattering with the %s, obtained via TRK scaling of rates A >= 12\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)

    # calc cumulative rate CDF for elastic scattering
    C1 = outer(zeros(len(gamma)),zeros(len(eps_elastic)))
    for i,x in enumerate(xs_elastic):
        C = iR.cumulative_rate_gamma_eps(eps_elastic,x,gamma,field)
        for i in range(0,len(C[:,0])):
            C[i,:] /= max(C[i,:])
        C1 += C
    C1 /= Niso
    C1[isnan(C1)] = 0.

    # save
    fname = 'data/ElasticScattering_CDF_%s.txt' % field.name
    output = c_[log10(gamma),C1]
    fmt = '%g' + '\t%g'*len(eps_elastic) 
    hdr = '# Elastic scattering with the %s averaged over all isotopes\n# log10(gamma), (1/lambda)_cumulative_normalized for eps = log10(2 keV) - log10(263 MeV) in 513 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)

    # Calculate branching ratios
    # for A < 12
    B1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1exc])
    for (Z, N, A) in isotopes1:
        s = (d1exc['Z'] == Z) * (d1exc['N'] == N)
        B1[s] /= sum(B1[s], axis=0)
    B1[isnan(B1)] = 0  # set to 0 when total cross section is 0

    # for A >= 12
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

    # Calculate photon emission probabilities
    B3 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs_exc_photon])
    PDchannel = transpose(array([d_sum_photon['Z'],d_sum_photon['N'],d_sum_photon['Zd'],d_sum_photon['Nd']]))
    for i,(Z, N, Zd, Nd) in enumerate(PDchannel):
        s = (d_exc_photon['Z'] == Z) * (d_exc_photon['N'] == N) * (d_exc_photon['Zd'] == Zd) * (d_exc_photon['Nd'] == Nd)
        B3[s] /= R3[i]
    B3[isnan(B3)] = 0  # set to 0 when total cross section (R3) is 0

    # save
    fname = 'data/pd_photon_emission_%s.txt'%field.name
    output = c_[d_exc_photon['Z'], d_exc_photon['N'], d_exc_photon['Zd'], d_exc_photon['Nd'], d_exc_photon['Ephoton']*1e6, B3]
    fmt = '%i\t%i\t%i\t%i\t%g' + '\t%g'*201
    hdr = 'Emission probabilities of photons with discrete energies via photo-disintegration with the %s\nZ, N, Z_daughter, N_daughter, Ephoton [eV], emission probability for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, output, fmt=fmt, header=hdr)
