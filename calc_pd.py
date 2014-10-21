from numpy import *
import photonField
import interactionRate as iR



eV = 1.60217657e-19
gamma = logspace(6, 14, 201)

# pad the cross-section data (extrapolate)
x  = logspace(5, 9, 501)[:-1] * eV  # [J]
eps = iR.romb_pad_logspaced(x, 1025)

dsum  = genfromtxt('pd-tables/xs_sumexcl_clean.txt')
dexcl = genfromtxt('pd-tables/xs_excl_thin.txt')
xs_sum  = array([iR.romb_pad_expo(x, y, eps) for y in  dsum[:,2:]]) * 1e-31  # [m^2]
xs_excl = array([iR.romb_pad_expo(x, y, eps) for y in dexcl[:,3:]]) * 1e-31  # [m^2]


fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Kneiske10(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Dole06(),
    photonField.EBL_Franceschini08(),
    photonField.CRB_Biermann96()
    ]


for field in fields:
    print field.name

    print '  calculate total interaction rate'
    R = array([iR.invMFP_fast(eps, y, gamma, field) for y in xs_sum])

    fname = 'data/pd_%s.txt' % field.name
    fmt = '%i\t%i' + '\t%g'*201
    hdr = 'Photo-disintegration with the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, c_[dsum[:,0:2], R], fmt=fmt, header=hdr)


    print '  calculate branching ratios'
    B = array([iR.invMFP_fast(eps, y, gamma, field) for y in xs_excl])

    # branching ratio = xs / sum(xs)
    for i in range(len(dsum)):
        Z, N = dsum[i,0:2]
        idx = (dexcl[:,0] == Z) * (dexcl[:,1] == N)
        B[idx] /= sum(B[idx], axis=0)

    B[isnan(B)] = 0  # set br to 1 when total cross section is 0

    fname = 'data/pd_branching_%s.txt'%field.name
    fmt = '%i\t%i\t%06d' + '\t%g'*201
    hdr = 'Photo-disintegration with the %s\nZ, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps' % field.info
    savetxt(fname, c_[dexcl[:,0:3], B], fmt=fmt, header=hdr)
