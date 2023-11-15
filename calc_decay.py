import numpy as np
from scipy.integrate import quad
import gitHelp as gh
import os
from units import mass_electron, c_squared, c_light, h_planck, keV, amu

cdir = os.path.split(__file__)[0]

# Script to preprocess the nuclear decay data table from the BNL NuDat2 database
# Decay Search: http://www.nndc.bnl.gov/nudat2/indx_sigma.jsp, output: formatted file --> decay_NuDat2.txt
# Decay Radiation Search: gamma_NuDat2.txt: http://www.nndc.bnl.gov/nudat2/indx_dec.jsp --> gamma_NuDat2.txt

class NuclearMassTable(object):
    """Class to provide tabulated nuclear masses
    
    This function reimplements the basic functionality of
    CRPropa's particleMass module see here:
    https://github.com/CRPropa/CRPropa3/blob/master/include/crpropa/ParticleMass.h
    """

    def __init__(self):
        try:
            datapath = os.path.join(cdir, 'data/nuclear_mass.txt')
            self.massTable = np.loadtxt(datapath, usecols=(2))
        except FileNotFoundError:
            print("The file 'data/nuclear_mass.txt' was not found.")
            print("Run the script calc_mass.py and try again.")

    def getMass(self, id: int) -> float:
        """Helper function to return tabulated nuclear masses
        
        id is not the usual CRPropa PID but id = Z * 31 + N
        with Z the charge number and N the neutron number.
        """
        return self.massTable[id]
    
    def nuclearMass(self, A: int, Z: int) -> float: 
        """nuclear mass for given mass (A) and charge (Z) number
        
        Particle masses that are not tabulated are approximated by
        A*amu-Z*mass_electron.
        """

        if ((A < 1) or (A > 56) or (Z < 0) or (Z > 26) or (Z > A)):
            print ("nuclearMass: nuclear mass not found in the mass table for A = {}, Z = {}. Approximated value used A * amu - Z * m_e instead.".format(A, Z))
            return A * amu - Z * mass_electron
        N = A - Z
        
        return self.getMass(Z * 31 + N)

class Decay:
    def load(self, s):
        """ extract decay parameter from a given line of the data file. """
        l = s.split('\t')
        self.Z = int(l[2])
        self.N = int(l[3])
        self.id = self.Z * 1000 + self.N

        # mode
        self.mode = l[12].strip()

        # decay time
        s = l[9].strip()
        if s == 'infinity':
            self.tau = np.inf
        elif s == '':
            self.tau = 0
        else:
            self.tau = float(s) / np.log(2)  # half-life --> life time

        # branching ratio
        s = ''.join(c for c in l[13] if c not in ('>','<','=','~','%',' ','?','\n'))
        self.brString = s
        if s == '':
            self.br = 0.
        else:
            self.br = float(s) / 100.  # % --> fraction

    def __str__(self):
        return 'Z=%i N=%i mode=%s tau=%.1e br=%.2f' % (self.Z, self.N, self.mode, self.tau, self.br)

    def isStable(self):
        """ returns if the nucleus is stable or not"""
        return self.tau == np.inf

    def isBetaPlus(self):
        """ returns if the nucleus has a beta plus decay mode"""
        return self.mode.find('E') > -1

    def isBetaMinus(self):
        """ returns if the nucleus has a beta minus decay mode"""
        return self.mode.find('B') > -1

class GammaEmission:
    def __init__(self, lines):
        l = lines[0].split('\t')
        self.Z = int(l[2])
        self.N = int(l[3])
        self.id = self.Z * 1000 + self.N
        self.mode = l[7].strip()

        self.energy = []
        self.intensity = []
        for line in lines:
            l = line.split('\t')
            self.energy.append(float(l[13]))
            self.intensity.append(float(l[17]))

    def __str__(self):
        s = 'Z = %i N = %i mode = %s' % (self.Z, self.N, self.mode)
        for i in range(len(self.energy)):
            s += '\n    energy = %.3f intensity = %.3e' % (self.energy[i], self.intensity[i])
        return s


### parse gamma emission data file
print ('\nParsing gamma emission data file')
print ('-------------------------------------')

datapath = os.path.join(cdir, 'tables/gamma_NuDat2.txt')
data = open(datapath)
lines = data.readlines()[1:-3] # skip header and footer
data.close()

# create list of gamma emission entries for each isotope
gammaTable = [[{} for n in range(31)] for z in range(27)]
for i, line in enumerate(lines):
    l = line.split('\t')
    Z = int(l[2])
    N = int(l[3])
    mode = l[7].strip()
    if (Z > 26) or (N > 30):  # skip if higher than Fe-56
        continue
    if (mode == 'IT'):  # skip isomeric transition
        continue
    if (l[4].strip() == '0+X' or float(l[4]) > 0):  # skip if parent nuclei in excited state
        continue
    if (l[11].strip() != 'G'):  # take only gamma radiation type
        continue
    if (l[12].strip() != ''):  # ionized nuclei: no Auger electrons, conversion electrons and annihilation
        continue
    gammaTable[Z][N].setdefault(mode, []).append(line)

# for each isotope and decay mode combine all gamma entries
for Z in range(27):
    for N in range(31):
        if not(gammaTable[Z][N]):  # no entry
            continue
        for mode, entries in gammaTable[Z][N].items():
            gammaTable[Z][N][mode] = GammaEmission(entries)


### explicitly edit some entries
print ('\nExplicitly editing certain entries')
print ('-------------------------------------')

# for beta-n decay of Na-27 photon emission probability is 100% if decay happens
g0 = gammaTable[11][16]['B-N']
g0.intensity[0] = 100.
print (g0, ' <- set photon emission probability to 100%\n')

# renormalize emission probability for beta+ decay of K-40 (BR = 10.86%, intensity = 10.66% -> emission prob if decay happens = 98.16%)
g0 = gammaTable[19][21]['EC']
g0.intensity[0] = 98.16
print (g0, ' <- renormalize photon emission probability to 98.16%\n')

# remove one of two tabulated beta- decays for K-46
g0 = gammaTable[19][27]['B-']
print (g0,'\n')
takeIndex = [2,3,4,7,9,11]
energy = []
intensity = []
for i in takeIndex:
    energy.append(g0.energy[i])
    intensity.append(g0.intensity[i])
g0.intensity = intensity
g0.energy = energy
print (g0, ' <- removed additional beta- decay with same properties\n')

# for beta- and beta+ decay of V-50 emission probability of photon is 100% if decay happens
g0 = gammaTable[23][27]['B-']
g1 = gammaTable[23][27]['EC']
g0.intensity[0] = 100.
g1.intensity[0] = 100.
print (g0, ' <- set photon emission probability to 100%\n')
print (g1, ' <- set photon emission probability to 100%\n')


### parse decay data file
print ('\nParsing decay data file')
print ('-------------------------------------')
datapath = os.path.join(cdir, 'tables/decay_NuDat2.txt')
fin = open(datapath)
lines = fin.readlines()
fin.close()

decayTable = [[[] for n in range(31)] for z in range(27)]

for line in lines[1:-3]:
    d = Decay()
    d.load(line)
    if (d.Z > 26) or (d.N > 30):
        continue
    if d.mode == 'IT':
        print (d, '<- skip (isomeric transition)')
        continue
    if d.tau == 0:
        print (d, '<- skip (missing lifetime)')
        continue
    if d.mode == '':
        if not(d.isStable()):
            print (d, '<- skip (missing decay mode)')
            continue
    decayTable[d.Z][d.N].append(d)


### remove duplicate decays
print ('\n\nRemoving duplicates')
print ('-------------------------------------')
for z in range(27):
    for n in range(31):
        dList = decayTable[z][n]

        if len(dList) < 2:
            continue

        for i, d1 in enumerate(dList):
            for d2 in dList[i+1:]:
                if d1.mode == d2.mode:
                    print (d1)
                    print (d2, ' <- remove \n')
                    dList.remove(d2)


### explicitly edit some entries
print ('\nExplicitly editing certain entries')
print ('-------------------------------------')

# remove Li-5 alpha decay (equivalent to existing proton emission)
d0 = decayTable[3][2][0]
d1 = decayTable[3][2][1]
print (d0)
print (d1, ' <- remove (equivalent to neutron emission)\n')
decayTable[3][2].remove(d1)

# remove He-5 alpha decay (equivalent to existing neutron emission)
d0 = decayTable[2][3][0]
d1 = decayTable[2][3][1]
print (d1)
print (d0, ' <- remove (equivalent to neutron emission)\n')
decayTable[2][3].remove(d0)

# modify B-12 "B3A" decay to "B2A" as it would leave an empty nucleus
d = decayTable[5][7][1]
print (d, ' <- change decay mode to B2A\n')
d.mode = 'B2A'

# Fe-45: to make beta+ decays exclusive
d = decayTable[26][19][0]
print (d, ' <- set branching ratio to 0 (ratio equal to sum of following ratios)')
d.br = 0
brSum = 0
for d in decayTable[26][19][1:]:
    print (d)
    brSum += d.br
for d in decayTable[26][19][1:]:
    d.br /= brSum


### calculate exclusive mean life times
print ('\n\nCalculating exclusive life times')
print ('-------------------------------------')
for z in range(27):
    for n in range(31):
        dList = decayTable[z][n]

        # skip for 0 or 1 entry
        if len(dList) < 2:
            continue

        # get sum of branching ratios
        brSum = 0
        for d in dList:
            brSum += d.br

        # if sum is 0, set branching ratios to equal values
        if brSum == 0:
            for d in dList:
                d.br = 1. / len(dList)

        # else if sum not 1, search for an inclusive decay and/or normalize the branching ratios
        elif brSum != 1.:
            dInclusive = None
            brSumExclusive = 0
            for i,d in enumerate(dList):
                if d.br == 1.0:
                    dInclusive = d # inclusive decay found
                else:
                    brSumExclusive += d.br # add exclusive branching ratio

            if dInclusive != None:
                if dInclusive.br <= brSumExclusive:
                    dList.remove(dInclusive) # remove if purely inclusive
                else:
                    dInclusive.br -= brSumExclusive # else make exclusive

            # normalize all branching ratios
            for d in dList:
                d.br /= brSum

        # finally, calculate exclusive decay time by dividing with branching ratio, while removing zeros
        for d in dList:
            if d.br == 0:
                print (d, ' <- remove (branching ratio 0)')
                dList.remove(d)
            else:
                d.tau /= d.br


### correct for electron capture contribution in beta+ decays
# see Basdevant, Fundamentals in Nuclear Physics, 4.3.2 and 4.3.3
print ('\nBeta+ correction')
print ('-------------------------------------')
Qe = mass_electron * c_squared  # electron energy [J]
a0 = 5.29177e-11  # Bohr radius [m]
hbar_c = c_light * (h_planck / 2 / np.pi)  # [m/J]

nucMass = NuclearMassTable()

for Z in range(27):
    for N in range(31):
        for d in decayTable[Z][N]:
            if not(d.isBetaPlus()):
                continue

            A  = Z+N
            m1 = nucMass.nuclearMass(A, Z)
            m2 = nucMass.nuclearMass(A, Z - 1)
            dm = (m1 - m2) * c_squared

            Qec   = (dm + Qe)
            Qbeta = (dm - Qe)

            # check if energetically possible
            if Qbeta < 0:
                print (d, ' <- make stable (beta+ decay not possible)')
                d.tau = np.inf
                continue

            f = lambda E: E * np.sqrt(E**2 - Qe**2) * (dm - E)**2
            I, Ierr = quad(f, Qe, dm)

            # ratio tau_beta+ / tau_ec
            f = np.pi**2 / 2 * (Z / a0*hbar_c)**3 * Qec**2 / I
            if f < 0:
                print (Qec)
            print (d, ' <- beta+ correction %.1e'%f)
            d.tau *= 1 + f

            # remove gamma decays which are not possible in beta+ decays
            try:
                g = gammaTable[Z][N]['EC']
            except:
                continue  # no gamma entry
            for i, Egamma in enumerate(g.energy):
                Egamma *= keV
                if Egamma > Qbeta:
                    print (d, ' <- remove gamma decay (Egamma %g < %g = Q)' % (Egamma, Qbeta))
                    g.energy.pop(i)
                    g.intensity.pop(i)
            if len(g.energy) == 0:
                gammaTable[Z][N].pop('EC')




### set immediate proton / neutron dripping for all other isotopes
print ('\n\nSet proton / neutron dripping for all other isotopes')
for z in range(0,27):
    for n in range(0,31):
        if (z + n)==0:
            continue

        dList = decayTable[z][n]

        if len(dList) > 0:
            continue

        # else set p/n dripping
        d = Decay()
        d.Z = z
        d.N = n
        d.tau = 1e-99
        d.br = 1.
        if z > n: # proton dripping
            d.mode = 'P'
        else: # neutron dripping
            d.mode = 'N'
        dList.append(d)

# output folder
folder = 'data'
if not os.path.exists(folder):
    os.makedirs(folder)
# save decay table
fout = open('data/nuclear_decay.txt','w')
# Add git hash of crpropa-data repository to header
try:
    git_hash = gh.get_git_revision_hash()
    fout.write('# Produced with crpropa-data version: '+git_hash+'\n')
except:
    pass 
fout.write('# Z, N, Decay Mode (#beta- #beta+ #alpha #p #n), Mean Life Time [s], Gamma Energy 1 [keV], Gamma Emission Probability 1, Gamma Energy 2 [keV], Gamma Emission Probability 2, ...\n')

# decay mode codes: #beta- #beta+ #alpha #p #n
modeDict = {'STABLE' : '0',
        'N'  : '00001',
        '2N' : '00002',
        'P'  : '00010',
        '2P' : '00020',
        'A'  : '00100',
        '2A' : '00200',
        'B-' : '10000',
        '2B-': '20000',
        'BN' : '10001',
        'B-N': '10001',
        'B2N': '10002',
        'B3N': '10003',
        'B4N': '10004',
        'BNA': '10101',
        'BA' : '10100',
        'B2A': '10200',
        'B3A': '10300',
        'B+' : '01000',
        'EC' : '01000',
        '2EC': '02000',
        'EA' : '01100',
        'EP' : '01010',
        'E2P': '01020',
        'E3P': '01030'}

for Z in range(27):
    for N in range(31):
        if (Z + N) == 0:
            continue

        for d in decayTable[Z][N]:
            if d.isStable():
                continue

            mode = modeDict[d.mode]
            s = '%i %i %s %e' % (Z, N, mode, d.tau)

            for key in gammaTable[Z][N].keys():
                if modeDict[key] != mode:
                    continue

                g = gammaTable[Z][N][key]
                for i in range(len(g.energy)):
                    s += ' %e %e'%(g.energy[i], g.intensity[i]/100.)

            fout.write(s + '\n')

fout.close()


### save isotopes with tau > 2 s to consider for photo-disintegration
# this is not needed for CRPropa
fout = open('data/isotopes-2s.txt', 'w')
# Add git hash of crpropa-data repository to header
try:
    git_hash = gh.get_git_revision_hash()
    fout.write('# Produced with crpropa-data version: '+git_hash+'\n')
except:
    pass 
fout.write('# Z\tN\tA\n')
fout.write('# isotopes with lifetime > 2s (including beta+ correction, see calc_decay.py)\n')
for z in range(1,27):
    for n in range(1,31):
        if (z + n)==0:
            continue
        c = 0  # total decay constant
        for d in decayTable[z][n]:
            c += 1 / d.tau
        if c < 1/2.:  # check for tau < 2 s
            fout.write('%i\t%i\t%i\n' % (z, n, z+n))
fout.close()


### save stable isotopes
# this is not needed for CRPropa
fout = open('data/isotopes-stable.txt', 'w')
# Add git hash of crpropa-data repository to header
try:
    git_hash = gh.get_git_revision_hash()
    fout.write('# Produced with crpropa-data version: '+git_hash+'\n')
except:
    pass 
fout.write('# Z\tN\tA\n')
fout.write('# stable isotopes (including beta+ correction, see calc_decay.py)\n')
for z in range(1,27):
    for n in range(1,31):
        if (z + n)==0:
            continue
        for d in decayTable[z][n]:
            if d.tau != np.inf:
                continue
            fout.write('%i\t%i\t%i\n' % (z, n, z+n))
fout.close()
