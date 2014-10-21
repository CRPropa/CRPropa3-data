from numpy import *
import os


def isBogus(Z, N, channel):
    s = channel  # expect string here
    dZ = int(s[1]) + int(s[2]) + int(s[3]) + 2*int(s[4]) + 2*int(s[5])
    dN = int(s[0]) + int(s[2]) + int(s[3]) +   int(s[4]) + 2*int(s[5])
    if (dZ + dN) == 0:
        print '    no photo-disintegration', channel
        return True
    if (dZ + dN) == (Z + N):
        print '    no nucleon left', channel
        return True
    if dZ > dZ:
        print '    too many protons lost', channel
        return True
    if dN > dN:
        print '    too many neutrons lost', channel
        return True
    return False


# -----------------------------------------------
# Collect all exclusive cross sections
# -----------------------------------------------
# output files
ftotal = open('xs_total.txt', 'w')
fsum   = open('xs_sumexcl.txt', 'w')
fexcl  = open('xs_excl.txt', 'w')

# output format
info = '#cross sections [mb] for incident photon energies eps = 10^5 - 10^8.992 eV\n'
fmt1 = '%i\t%i\t%s' + '\t%g'*500 + '\n'
fmt2 = '%i\t%i' + '\t%g'*500 + '\n'

ftotal.write('# Z\tN\txs\n' + info)
fsum.write(  '# Z\tN\txs\n' + info)
fexcl.write( '# Z\tN\tchannel\txs\n' + info)

dxs = {}
base = 'data/'

for folder in os.listdir(base):
    print folder
    s = folder.split('-')
    Z, A = int(s[0]), int(s[1])
    N = A - Z

    xs = genfromtxt(base + folder + '/total.tot', usecols=1)
    ftotal.write(fmt2 % ((Z, N) + tuple(xs)))

    # loop over exclusive channels and add to dictionary
    dxs[Z, N] = []
    for f in os.listdir(base + folder):
        if not( f.startswith('xs') and f.endswith('tot') ):
            continue  # consider the xs....tot files only

        channel = f.strip('xs.tot')
        if isBogus(Z, N, channel):
            continue  # skip bogus channels

        fname = os.path.join(base, folder, f)
        xs = genfromtxt(fname, usecols=1)
        dxs[Z, N].append((channel, xs))


# -----------------------------------------------
# Save exclusive and sum of exclusive cross sections sorted by element
# -----------------------------------------------
keys = dxs.keys()
keys.sort()
for k in keys:
    xstot = zeros(500)
    Z, N = k
    for ch_xs in dxs[k]:
        ch, xs = ch_xs
        fexcl.write( fmt1 % ((Z, N, ch) + tuple(xs)) )
        xstot += xs  # sum up exclusive channels
    fsum.write( fmt2 % ((Z, N) + tuple(xstot)) )

ftotal.close()
fsum.close()
fexcl.close()
