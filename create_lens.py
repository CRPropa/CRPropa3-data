#!/usr/bin/env python
from __future__ import print_function, division
import crpropa
import h5py
import healpy
import scipy.sparse
import numpy as np
import os
import sys
import time
import struct
import logging

import multiprocessing
import argparse
import gitHelp as gh

log = logging.getLogger("create_lens")
ch = logging.StreamHandler()
log.addHandler(ch)

parser = argparse.ArgumentParser(description="""
Creates a lens from backtracking simulations of anti-particles in the galactic
magnetic field. See the CRPropa examples for an setup for such backtracking
simulations. Simulation files are assumed to contain only anti-protons, every
file may contain arbitrary energies. Lens creation accounts for any
anisotropies in particle emission e.g. from random starting directions. Please
add additional information on how the lens was created to a README file in the
lens directory.
""", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--rmin", default=17.00, help="minimum rigidity of lens log10(Rmin / eV)")
parser.add_argument("--rmax", default=21.00, help="maximum rigidity of lens log10(Rmax / eV)")
parser.add_argument("--dr", default=0.1, help="width of rigidity bins (10**Ri eV... 10**Ri+dr eV)")
parser.add_argument("--nproc", default=multiprocessing.cpu_count(), help="number of processes. Default is one per available CPU.")
parser.add_argument("-v", action="store_true", help="Print debug output")
parser.add_argument("lensname", help="name of the lens / Output directory the lens is created in.")
parser.add_argument("inputfiles", help="input files to create lens from", nargs="+")

args = parser.parse_args()

if args.v:
    print("Set log level to DEBUG")
    log.setLevel(logging.DEBUG)

if os.path.exists(args.lensname):
    print("Error - path '{}' already exists!".format(args.lensname))
    exit(-1)


print('Creating lens {} from {} input files using {} processes:'.format(args.lensname, len(args.inputfiles), args.nproc))
os.mkdir(args.lensname)

rigidities = np.arange(args.rmin, args.rmax, args.dr)
# Hardcoded here as also in the lenses
nside = 64
npix = healpy.nside2npix(nside)

matrices = {}

tot_cr = 0
print(' - Reading data:', end='')
if args.v:
    print()
else:
    sys.stdout.flush()


start = time.time()
for i, filename in enumerate(args.inputfiles):
    log.debug("  Open File {}".format(filename))
    infile = h5py.File(filename)
    data = infile['CRPROPA3'][:]
    log.debug("  Data size: {}".format(data.size))

    log.debug("  Calculating out pixel ...")
    pix_out = healpy.vec2pix(nside, data['Px'], data['Py'], data['Pz'])
    log.debug("  Calculating in pixel ...")
    pix_in  = healpy.vec2pix(nside, data['P0x'], data['P0y'], data['P0z'])
    tot_cr += pix_in.size
    log.debug("  Associate data to rigidity matrices ...")
    for rmin in rigidities:
        idx = (data['E']*1E18 >= 10**rmin) * (data['E']*1E18 < 10**(rmin+args.dr))
        log.debug("   lg R / eV = {} : {} cosmic rays ".format(rmin, idx.sum()))
        M = scipy.sparse.coo_matrix((np.ones_like(pix_in[idx]), (pix_in[idx], pix_out[idx])), shape= [npix, npix])
        if rmin not in matrices:
            matrices[rmin] = scipy.sparse.coo_matrix((npix, npix))
        matrices[rmin] += M


    print("\r - Reading data: read {:3.0f}%, got {} cosmic rays in total ({:.1f} per rigidity and pixel).".format( (i + 1) / len(args.inputfiles) * 100, tot_cr, tot_cr / (npix * rigidities.size )), end='')
    if args.v:
        print()
    else:
        sys.stdout.flush()

    infile.close()
print()

duration = time.time() - start
log.debug("Total time to read data: {}s, {}s per file".format(duration, duration / len(args.inputfiles)))


cfg = open(os.path.join(args.lensname, 'lens.cfg'), 'w')
cfg.write('# Magnetic lens {} for usage with CRPRopa3\n'.format(args.lensname))
cfg.write('# Created on {} with {} backtracked particles\n'.format(time.ctime(), tot_cr))
try:
    cfg.write('# Produced with crpropa-data version: {}\n'.format(gh.get_git_revision_hash())
except:
    pass
cfg.write('#\n')
cfg.write('# fname log(Rmin / eV) log(Rmax / eV)\n')
cfg.write('#\n')

def rigidity_processor(filename, data):
    try:
        mldat_of = open(filename, 'wb')
        # Write header information
        mldat_of.write(struct.pack("<I", data.nnz))
        mldat_of.write(struct.pack("<I", npix))
        mldat_of.write(struct.pack("<I", npix))

        # Normalize rows to account for anisotropic particle emission
        K = data.tocsr()
        for i, row in enumerate(K):
            S = row.sum()
            if S > 0:
                K[i] /= S

        # Write output
        M = K.tocoo()
        Iv, Jv = M.nonzero()
        for i, d in enumerate(M.data):
            mldat_of.write(struct.pack("<I", Iv[i]))
            mldat_of.write(struct.pack("<I", Jv[i]))
            mldat_of.write(struct.pack("<d", d))
        mldat_of.close()
    except Exception as E:
        logging.exception(E)
        return False
    return True

pool = multiprocessing.Pool(processes=args.nproc)

results = []
for rmin in rigidities:
    ofname =  "{:.2f}.mldat".format(rmin)
    cfg.write('{:20} {:.2f}     {:.2f}\n'.format(ofname, rmin, rmin+args.dr))
    filename = os.path.join(args.lensname, ofname)
    pool.apply_async(rigidity_processor, ((filename, matrices[rmin])), callback=results.append)
cfg.close()

pool.close()
last_len = -1
while len(results) != len(rigidities):
    if len(results) != last_len:
        last_len = len(results)
        print("\r - Normalizing data and writing output files: {:3.0f}% finished".format(len(results) / len(rigidities) * 100 ), end='')
        sys.stdout.flush()
    time.sleep(0.5)
print("\r - Normalizing data and writing output files: {:3.0f}% finished".format(len(results) / len(rigidities) * 100 ))
duration = time.time() - start
log.debug("Total time needed {}s.".format(duration))



