import numpy as np
import gitHelp as gh
import os

from units import amu, mass_electron, mass_proton, mass_neutron

cdir = os.path.split(__file__)[0]

# This script generates a table of nuclear mass for all combinations (Z,N) Z=0..26, N=0..30
# For measured atoms, the NIST data table is used.
# All other atomic masses are taken to be A * amu
# The nuclear mass is then approximated as atomic mass minus electron mass.
# Here, electron binding energies (~keV) can be neglected compared to nuclear binding energies (~MeV) and nucleon and electron masses (~GeV, ~MeV).

### read NIST data
# See: http://www.nist.gov/pml/data/comp.cfm
# All Isotopes, Linearized ASCII Output
# http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all
datapath = os.path.join(cdir, 'tables/mass_NIST.txt')
fin = open(datapath, 'r')
D = np.zeros((27, 31))

for i in range(4):
    fin.readline() # skip header

for line in fin.readlines():
    if line.startswith('Atomic Number'):
        z = int(line.strip('Atomic Number = '))
        continue

    if line.startswith('Mass Number'):
        a = int(line.strip('Mass Number = '))
        continue

    if line.startswith('Relative Atomic Mass'):
        line = line.strip('Relative Atomic Mass = ')
        relAtomicMass = line[:line.find('(')]
        n = a - z
        if a == 1:
            continue # skip H-1
        if z > 26:
            continue # skip isotopes with Z > 26
        if n > 30:
            continue # skip isotopes with N > 30

        # mass in [kg] minus mass of electrons
        D[z, n] = float(relAtomicMass) * amu - z * mass_electron


### add neutron and proton mass
D[1, 0] = mass_proton
D[0, 1] = mass_neutron

### fill empty entries in table with A * amu - Z * m_e approximation
for z in range(27):
    for n in range(31):
        if D[z, n] == 0:
            D[z, n] = (z + n) * amu - z * mass_electron


# output folder
folder = 'data'
if not os.path.exists(folder):
    os.makedirs(folder)
# Write to file
fout = open('data/nuclear_mass.txt', 'w')

# Add git hash of crpropa-data repository to header
try:
    git_hash = gh.get_git_revision_hash()
    fout.write('# Produced with crpropa-data version: '+git_hash+'\n')
except:
    pass 

fout.write('# Nuclear Mass of Isotopes Z < 26 and A < 56\n')
for z in range(27):
    for n in range(31):
        fout.write(str(z) + ' ' + str(n) + ' ' + str(D[z, n]) + '\n')

fout.close()