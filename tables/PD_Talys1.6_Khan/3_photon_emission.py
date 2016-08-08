from numpy import *
import os

# -----------------------------------------------
# Collect exclusive gamma cross sections from TALYS for remaining files after thinning
# -----------------------------------------------

fexcl = open('xs_excl_photon.txt', 'w')
fexcl.write('# Z\tN\tZ_daughter\tN_daughter\tEgamma [MeV]\txs [mb]\n')
fexcl.write('#cross sections [mb] for emission of photon with Egamma [MeV] for incident background photon energies eps = 0.2 - 200 MeV in steps of logE = 0.01\n')
fmt = '%i\t%i\t%i\t%i\t%.4f' + '\t%.4g'*301 + '\n'  # output format

fthin = genfromtxt('xs_thin.txt')
n = len(fthin[:,0])
thin = outer(zeros(n),zeros(4))
dat = open('xs_thin.txt')
for i,line in enumerate(dat):
    l = line.split()
    Z = int(l[0])
    N = int(l[1])
    channel = l[2]
    Zd = Z - (int(channel[1]) + int(channel[2]) + int(channel[3]) + 2*int(channel[4]) + 2*int(channel[5]))
    Nd = N - (int(channel[0]) + int(channel[2]) + int(channel[3]) + int(channel[4]) + 2*int(channel[5]))
    thin[i,0] = Z 
    thin[i,1] = N
    thin[i,2] = Zd
    thin[i,3] = Nd

dat.close()
isotopes_part1 = genfromtxt('../PD_external/isotopes.txt') # note: TALYS can not process H, He
isotopes_part2 = genfromtxt('isotopes.txt')
isotopes = vstack((isotopes_part1, isotopes_part2))

for (Z,N,A) in isotopes:
    print Z, N
    folder = '%i-%i/' % (Z, A)

    for f in os.listdir(folder):
        if not( f.startswith('gam') and f.endswith('tot') ):
            continue  # consider the gam....tot files only

        channel = f.strip('gam.tot')
        Zd = int(channel[:3])
        Nd = int(channel[3:6]) - Zd
        check = (Z == thin[:,0]) * (N == thin[:,1]) * (Zd == thin[:,2]) * (Nd == thin[:,3]) # check if PD channel remains after xs file thinning
        if sum(check) == 0:
            continue
        Egamma = open(folder + f).readlines()[0].split()
        Egamma = double(Egamma[len(Egamma)-1])
        xs = genfromtxt(folder + f, usecols=1)
        fexcl.write( fmt % ((Z, N, Zd, Nd, Egamma) + tuple(xs)) )

fexcl.close()

# get total cross section for all PDs from mother isotope to daughter isotope
fdaughter = open('xs_daughter.txt', 'w')
fmt2 = '%i\t%i\t%i\t%i' + '\t%.4g'*301 + '\n'  # output format

for (Z,N,A) in isotopes:
    print Z, N
    folder = '%i-%i/' % (Z, A)

    for f in os.listdir(folder):
        if not( f.startswith('xs') and f.endswith('tot') ):
            continue  # consider the xs....tot files only

        channel = f.strip('xs.tot')
        Zd = Z - (int(channel[1]) + int(channel[2]) + int(channel[3]) + 2*int(channel[4]) + 2*int(channel[5]))
        Nd = N - (int(channel[0]) + int(channel[2]) + int(channel[3]) + int(channel[4]) + 2*int(channel[5]))
        check = (Z == thin[:,0]) * (N == thin[:,1]) * (Zd == thin[:,2]) * (Nd == thin[:,3]) # check if PD channel remains after xs file thinning
        if sum(check) == 0:
            continue

        xs = genfromtxt(folder + f, usecols=1)
        fdaughter.write( fmt2 % ((Z, N, Zd, Nd) + tuple(xs)) )

fdaughter.close()

dat = genfromtxt('xs_daughter.txt')
fsum_daughter = open('xs_sum_daughter.txt', 'w')
fsum_daughter.write('# Z\tN\tZ_daughter\tN_daughter\txs [mb]\n')
fsum_daughter.write('#total cross sections [mb] for all PD from mother isotope to daughter isotope for photon energies eps = 0.2 - 200 MeV in steps of logE = 0.01\n')
PDs = vstack({tuple(r) for r in dat[:,:4]})
PDs = PDs.astype(int)
PDs = sort(PDs.view('i8,i8,i8,i8'), order=['f0','f1','f2','f3'], axis=0).view(int)

for (Z, N, Zd, Nd) in PDs:
    idx = (Z == dat[:,0]) * (N == dat[:,1]) * (Zd == dat[:,2]) * (Nd == dat[:,3])
    xs = dat[:,4:][idx]
    xs_sum = sum(xs, axis=0)
    fsum_daughter.write( fmt2 % ((Z, N, Zd, Nd) + tuple(xs_sum)) )

fsum_daughter.close()
os.remove('xs_daughter.txt')
