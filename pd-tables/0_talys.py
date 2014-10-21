from numpy import *
import os
import time

elements = {
    1  :  'H',
    2  : 'He',
    3  : 'Li',
    4  : 'Be',
    5  :  'B',
    6  :  'C',
    7  :  'N',
    8  :  'O',
    9  :  'F',
    10 : 'Ne',
    11 : 'Na',
    12 : 'Mg',
    13 : 'Al',
    14 : 'Si',
    15 :  'P',
    16 :  'S',
    17 : 'Cl',
    18 : 'Ar',
    19 :  'K',
    20 : 'Ca',
    21 : 'Sc',
    22 : 'Ti',
    23 :  'V',
    24 : 'Cr',
    25 : 'Mn',
    26 : 'Fe'}

input_str = """projectile g
energy eps
element %s
mass %i
channels y
maxchannel 8
isomer 1.e38
fileresidual n
components n
"""

eps = logspace(-1, 3, 501)[:-1]  # MeV

d = genfromtxt('isotopes.txt')
for z, a in d[94:]:
    print z, a
    time.sleep(3)

    folder = '%i-%i' % (z, a)
    try:
        os.mkdir(folder)
    except:
        pass

    os.chdir(folder)

    f = open('input', 'w')
    f.write(input_str % (elements[z], a))
    f.close()
    savetxt('eps', eps, fmt='%.6g')

    os.system('/user/walz/local/bin/talys < input > output')
    os.chdir('..')
