from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19
gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors

fields = [
    photonField.CMB(),
    photonField.EBL_Gilmore12()]


# Load cross section data from TALYS
ddir = 'tables/PD_Talys1.8_Khan/'
eps  = genfromtxt(ddir+'eps_elastic.txt') * eV * 1e6  # nuclear rest frame photon energies [J]
data = genfromtxt(ddir+'xs_elastic.txt', dtype=[('Z',int), ('N',int), ('xs','%if8'%len(eps))])

# Only consider TALYS cross sections for A >= 12
idx = (data['Z'] + data['N']) >= 12
data = data[idx]

# Factor out the principal scaling given by the TRK formula: sigma_int ~ Z*N/A
data['xs'] /= ( data['Z'] * data['N'] / (data['Z'] + data['N']) )[:,newaxis]

# Pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration
eps = iR.romb_pad_logspaced(eps, 513)
xs  = array([iR.romb_pad_zero(x, 513) for x in data['xs']]) * 1e-31



for field in fields:
    print (field.name)

    # output folder
    folder = 'data/ElasticScattering'
    if not os.path.exists(folder):
        os.makedirs(folder)

    # Calculate the interaction rate, averaged over all isotopes
    rate = mean([iR.calc_rate_eps(eps, x, gamma, field) for x in xs], axis=0)

    savetxt(folder + '/rate_%s.txt' % field.name.split('_')[0],
        rate, fmt='%g',
        header='Average interaction rate for elastic scattering of %s photons off nuclei\nScale with Z*N/A for nuclei\n1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info)

    # Calculate CDF for background photon energies, averaged over all isotopes
    CDF = zeros((len(gamma), len(eps)))
    for x in xs:
        C = iR.calc_diffrate_eps(eps, x, gamma, field)
        CDF += C / amax(C, axis=1)[:,newaxis]
    CDF /= len(data)
    CDF = nan_to_num(CDF)

    savetxt(folder + '/cdf_%s.txt' % field.name.split('_')[0],
        c_[log10(gamma), CDF],
        fmt='%g' + '\t%g'*len(eps),
        header='# Average CDF(background photon energy) for elastic scattering with the %s\n# log10(gamma), (1/lambda)_cumulative for eps = log10(2 keV) - log10(263 MeV) in 513 steps' % field.info)
