from numpy import *

# -------------------------------------------------
# Thin out photo-disintegration channels: For each isotope, select only those
# channels, whose cross section is at least 1% of the total cross section
# for any photon energy
# -------------------------------------------------

d1 = genfromtxt('xs_sumexcl.txt')
d2 = genfromtxt('xs_excl.txt')

select = ones(len(d2), dtype=bool)  # mask of selected channels
for i in range(len(d1)):
    Z, N, xstot = d1[i,0], d1[i,1], d1[i,2:]
    idx = (d2[:,0] == Z) * (d2[:,1] == N)
    xs = d2[idx,3:]
    ratio = xs / xstot  # relativ cross section
    ratio[isnan(ratio)] = 0  # remove nans
    select[idx==True] *= amax(ratio, axis=1) > 0.01  # store selection

print 'selected', sum(select), 'of' ,len(d2), 'channels'

fmt = '%i\t%i\t%06d\t' + '%g\t'*500
savetxt('xs_excl_thin.txt', d2[select], fmt=fmt)
