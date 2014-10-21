from pylab import *

# -------------------------------------------------
# Remove the constant offsets at lower energies, seen for a few isotopes
# -------------------------------------------------

dt = [('Z','u1'), ('N','u1'), ('xs','500f8')]
d = genfromtxt('xs_sumexcl.txt', dtype=dt)
d.sort()

for xs in d['xs']:
    # first non-zero element
    i0 = nonzero(xs)[0][0]
    # last occurrence of the first non-zero element
    i1 = searchsorted(xs[i0:], xs[i0], side='right')
    if i1 - i0 > 5:  # remove if constant plateau
        xs[i0:i1+1] = 0

with open('xs_sumexcl_clean.txt','w') as f:
    fmt = '%i\t%i' + '\t%g'*500 + '\n'
    for Z, N, xs in d:
        f.write( fmt % ((Z, N) + tuple(xs)) )
