from numpy import *
import os

# -----------------------------------------------
# Collect exclusive photon emission cross sections from TALYS for remaining files after thinning
# -----------------------------------------------

data = genfromtxt('xs_excl_photon.txt')
s = ones(len(data[:,0]), dtype=bool)
daughter = vstack({tuple(r) for r in data[:,:4]})
daughter = daughter.astype(int)
daughter = sort(daughter.view('i8,i8,i8,i8'), order=['f0','f1','f2','f3'], axis=0).view(int)

for (Z,N,Zd,Nd) in daughter:
    print Z, N, Zd, Nd
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    xs = data[:,5:][idx]
    xs_sum = sum(xs, axis=0)
    # calculate branching ratios, set to 0 when total cross-section = 0
    br = xs / xs_sum
    br[isnan(br)] = 0

    # find channels with branching ratios above the threshold
    s[idx] = amax(br, axis=1) > 0.05

print 'selected', sum(s), 'of' ,len(data), 'channels'
savetxt('xs_photon_thin.txt', data[s], fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
