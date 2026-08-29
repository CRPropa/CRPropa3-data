#!/usr/bin/env python3
"""Export the AstroPhoMes empirical photomeson model to CRPropa data tables.

The durable specification of these files is DATA_GENERATION.txt in this
repository.

Reads the AstroPhoMes checkout (env ASTROPHOMES_PATH, or --astrophomes) and
writes three field-independent tables:

    basis.txt           159-point eps_r grid and the cross-section basis vectors
    redistribution.txt  SOPHIA x-spectra (truncated, cumulative)
    fragments.txt       exclusive fragmentation channels

Empirical model: Morejon, Fedynitch, Boncioli, Biehl & Winter,
JCAP 11 (2019) 007, arXiv:1904.07999.

Usage:
    python3 export_tables.py --out <build>/data/PhotoPionProductionEmpirical
"""

import argparse
import importlib.util
import os
import sys

import numpy as np

# products kept from the redistribution tables.  The muon (7, 10) and neutrino
# (11..16) ids are dropped: the tables carry UNDECAYED mesons, so those entries
# are identically zero (verified) and the decay chain is done in the module.
PRODUCTS = [
    (0, "gamma"), (2, "pi+"), (3, "pi-"), (4, "pi0"),
    (20, "e-"), (21, "e+"), (50, "K+"), (51, "K-"),
    (100, "n"), (101, "p"),
]

# light fragment basis, in the order used by fragments.txt.  Deliberately
# excludes AstroPhoMes' unphysical species (202 diproton, 303/403 Li-4).
FRAGMENTS = [(402, 4, 2), (302, 3, 2), (301, 3, 1), (201, 2, 1), (101, 1, 1), (100, 1, 0)]

# Drop residual channels below this weight.  CRISP uses 1e-3, which is fine when
# the channels are folded into rates, but here <dA> is the observable and that cut
# costs 2.9% of it (5.649 vs 5.815 for Fe-56) by discarding the deep-spallation
# tail.  1e-4 recovers <dA> to 5.805 for ~40% more rows.
MULT_CUT = 1e-4
XSEC_CUT = 1e-8   # drop redistribution rows below this multiplicity

# Default path to the AstroPhoMes submodule, relative to this file's location
_DEFAULT_ASTROPHOMES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "../../lib/AstroPhoMes")


def load_astrophomes(repo):
    """Import AstroPhoMes, whose modules do `from config import *`."""
    repo = os.path.abspath(repo)
    if not os.path.isdir(repo):
        raise SystemExit("AstroPhoMes not found at %s (set ASTROPHOMES_PATH)" % repo)
    spec = importlib.util.spec_from_file_location("config", os.path.join(repo, "config.py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules["config"] = module
    spec.loader.exec_module(module)
    if repo not in sys.path:
        sys.path.insert(0, repo)
    models = importlib.import_module("photomeson_lib.photomeson_models")
    phenom = importlib.import_module("photomeson_lib.phenom_relations")
    return models, phenom


def write_basis(em, path):
    """eps_r grid plus the vectors from which sigma_nonel and G are rebuilt.

    All cross sections are in microbarn, including pion_spl, which is compared
    against the nucleon cross sections in AstroPhoMes' cs_incl_diff (it is
    evaluated at eps_r in MeV, as AstroPhoMes does).  alpha and alpha_pi are the
    mass-scaling exponents of the nonelastic and the pion cross section.
    """
    e = em.egrid
    cols = np.column_stack([
        e,
        1e30 * em.univ_spl(e),      # cm^2 -> microbarn
        em.cs_proton_grid,
        em.cs_neutron_grid,
        em.alpha(e),
        em.alphapi_spl(e),
        em.pion_spl(e * 1e3),       # spline is parameterised in MeV
    ])
    with open(path, "w") as f:
        f.write("# Empirical photomeson model basis vectors (Morejon et al. 2019, arXiv:1904.07999)\n")
        f.write("# Produced by dataGeneration/PhotoPionProductionEmpirical/export_tables.py\n")
        f.write("#\n")
        f.write("# WARNING: the 159-point eps_r grid is load-bearing and must not be re-gridded.\n")
        f.write("# AstroPhoMes' pion fade windows are hardcoded to the INDICES 0..31 and 55..104\n")
        f.write("# of this grid (photomeson_models.py:358-359, \"hardcoded, found manually\").\n")
        f.write("# The grid is not log-equidistant: dlog10 ranges from 0.0092 to 0.197.\n")
        f.write("#\n")
        f.write("# pion_spl is an extrapolating spline and goes NEGATIVE outside its fit range\n")
        f.write("# (203-788 MeV); it is rescued only by vanishing fade weights.  Consumers must\n")
        f.write("# clip the blended pion scaling G to >= 0.\n")
        f.write("#\n")
        f.write("# eps_r[GeV]  sigma_univ[ub]  sigma_p[ub]  sigma_n[ub]  alpha  alpha_pi  pion_spl[mb]\n")
        for row in cols:
            f.write("%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\n" % tuple(row))
    return len(e)


def write_redistribution(em, path, renormalise=True):
    """SOPHIA x-spectra, truncated to the contiguous nonzero range, as CDFs.

    Only 26.6% of cells are nonzero and the nonzero range at fixed
    (product, eps_r) is contiguous, so storing [first, first+n) shrinks the file
    ~2.6x.  Each row carries the bin-summed multiplicity and a normalised
    cumulative, which is what the module samples by inverse-CDF.

    With renormalise=True the multiplicities at each eps_r are scaled by a
    single common factor so that sum_s int x dN_s/dx dx == 1, i.e. the kept
    species carry exactly one struck nucleon's energy.  The factor is recorded
    per eps_r bin in the header.  A common factor leaves the p/n charge ratio
    (the only way nucleon multiplicities are used) untouched.
    """
    e = em.egrid
    xb = em.xbins
    xw = xb[1:] - xb[:-1]
    xc = 0.5 * (xb[1:] + xb[:-1])
    nEps, nX = len(e), len(xw)
    redist = {1: em.redist_proton, 0: em.redist_neutron}

    # energy carried by the kept species, per nucleon and eps_r bin
    carried = np.zeros((2, nEps))
    for nucleon in (0, 1):
        for pid, _ in PRODUCTS:
            carried[nucleon] += (redist[nucleon][pid] * xw * xc).sum(axis=1)
    scale = np.ones_like(carried)
    if renormalise:
        scale = np.where(carried > 0, 1.0 / np.where(carried > 0, carried, 1.0), 1.0)

    rows = []
    for nucleon in (0, 1):
        for pid, _ in PRODUCTS:
            tab = redist[nucleon][pid]
            for k in range(nEps):
                dndx = tab[k]
                w = dndx * xw
                mult = w.sum() * scale[nucleon, k]
                if mult < XSEC_CUT:
                    continue
                nz = np.nonzero(w)[0]
                first, last = int(nz[0]), int(nz[-1])
                seg = w[first:last + 1]
                cdf = np.cumsum(seg) / seg.sum()
                cdf[-1] = 1.0
                rows.append((nucleon, pid, k, first, len(seg), mult, cdf))

    with open(path, "w") as f:
        f.write("# SOPHIA redistribution tables x = E_secondary / E_nucleon, from AstroPhoMes\n")
        f.write("# Produced by dataGeneration/PhotoPionProductionEmpirical/export_tables.py\n")
        f.write("#\n")
        f.write("# These spectra are UNDECAYED: pi+-, pi0 and K+- appear as such, the neutrino\n")
        f.write("# entries are identically zero and the gamma/e+- entries are small direct\n")
        f.write("# emission channels.  The decay chain is the consumer's responsibility.\n")
        f.write("#\n")
        f.write("# x bins: %d edges, log-equidistant, log10(x) = %.15f .. %.15f\n"
                % (len(xb), np.log10(xb[0]), np.log10(xb[-1])))
        f.write("# eps_r grid: the %d points of basis.txt, referenced by index\n" % nEps)
        f.write("# products: %s\n" % ", ".join("%d=%s" % p for p in PRODUCTS))
        f.write("# energy renormalisation: %s\n" % ("on" if renormalise else "off"))
        f.write("# sum_s <x_s> before renormalisation, proton parent, per eps_r index:\n")
        for k in range(0, nEps, 20):
            f.write("#   i=%3d eps_r=%.4e  %.6f\n" % (k, e[k], carried[1, k]))
        f.write("#\n")
        f.write("# nucleon(1=p,0=n)  product  epsIndex  first  n  multiplicity  cdf[n]\n")
        for nucleon, pid, k, first, n, mult, cdf in rows:
            f.write("%d\t%d\t%d\t%d\t%d\t%.9e\t%s\n"
                    % (nucleon, pid, k, first, n, mult, "\t".join("%.7e" % v for v in cdf)))
    return len(rows)


def decompose(dA, dZ, mix):
    """Split (dA, dZ) into counts of [He4, He3, t, d, p, n].

    Clusters are allocated in proportion to the mother's own light-fragment mix
    (scaled so the cluster mass sums to dA), then free protons and neutrons close
    A and Z exactly.  Clusters are backed off whenever they would overdraw the
    available charge or neutron number, so the result always closes.
    """
    dN = dA - dZ
    if dA <= 0 or dZ < 0 or dN < 0:
        return None
    comp = [(a, z) for (_, a, z) in FRAGMENTS[:4]]          # He4, He3, t, d
    mass = sum(m * a for m, (a, _) in zip(mix[:4], comp))
    counts = [0, 0, 0, 0]
    if mass > 0:
        target = [m * dA / mass for m in mix[:4]]
        order = sorted(range(4), key=lambda i: -target[i])
        useA, useZ, useN = 0, 0, 0
        for i in order:
            a, z = comp[i]
            n = int(round(target[i]))
            while n > 0:
                if (useA + n * a <= dA and useZ + n * z <= dZ
                        and useN + n * (a - z) <= dN):
                    counts[i] = n
                    useA += n * a
                    useZ += n * z
                    useN += n * (a - z)
                    break
                n -= 1
    usedZ = sum(c * z for c, (_, z) in zip(counts, comp))
    usedN = sum(c * (a - z) for c, (a, z) in zip(counts, comp))
    return counts + [dZ - usedZ, dN - usedN]


def write_fragments(models, phenom, em, path):
    """Exclusive fragmentation channels, one heavy residual per event.

    The empirical multiplicity table is INCLUSIVE (total multiplicity ~4.2 for
    Fe-56).  Following CRISP's load_astrophomes(channels='empirical'), keep only
    the heavy survivors, drop negligible cells and renormalise to exactly one
    event, then decompose each channel's (dA, dZ) into light clusters.  These
    weights are energy-independent (the empirical formulas are energy averaged)
     so a single scalar per (mother, residual) suffices.
    """
    mothers = sorted(k for k in models.spec_data
                     if isinstance(k, int) and k >= 200)
    rows, stats = [], {}
    for mom in mothers:
        A, Z = mom // 100, mom % 100
        N = A - Z
        if Z < 1 or N < 0:
            continue
        try:
            mults = phenom.multiplicity_table(mom)
        except Exception as exc:                      # noqa: BLE001 - report and skip
            print("  skipping %d: %s" % (mom, exc))
            continue

        mix = [float(mults.get(fid, 0.0)) for fid, _, _ in FRAGMENTS]
        lo = max(1, (A + 1) // 2)
        keep = {}
        for dau, w in mults.items():
            Ad, Zd = dau // 100, dau % 100
            if not (A > Ad >= lo) or Zd < 0 or Ad - Zd < 0:
                continue
            w = float(w)
            if w >= MULT_CUT:
                keep[(Zd, Ad - Zd)] = keep.get((Zd, Ad - Zd), 0.0) + w
        if not keep:
            continue
        total = sum(keep.values())

        dAsum = 0.0
        for (Zd, Nd), w in sorted(keep.items()):
            w /= total
            dec = decompose(A - (Zd + Nd), Z - Zd, mix)
            if dec is None:
                continue
            rows.append((Z, N, Zd, Nd, dec, w))
            dAsum += (A - (Zd + Nd)) * w
        stats[mom] = (len(keep), total, dAsum)

    with open(path, "w") as f:
        f.write("# Empirical photomeson model exclusive fragmentation channels\n")
        f.write("# Produced by dataGeneration/PhotoPionProductionEmpirical/export_tables.py\n")
        f.write("#\n")
        f.write("# One heavy residual per event (A_residual >= max(1, ceil(A/2))); cells below\n")
        f.write("# %g of an event are dropped and the remainder renormalised to sum 1 per\n" % MULT_CUT)
        f.write("# mother.  Weights are energy independent.  Light fragments close A and Z\n")
        f.write("# exactly, so mass and charge are conserved per event by construction.\n")
        f.write("#\n")
        if 5626 in stats:
            n, raw, dA = stats[5626]
            f.write("# Fe-56: %d residual channels, raw sum %.4f, <dA> = %.4f\n" % (n, raw, dA))
        f.write("#\n")
        f.write("# Z  N  Zd  Nd  nHe4 nHe3 nH3 nH2 nP nN  weight\n")
        for Z, N, Zd, Nd, dec, w in rows:
            f.write("%d\t%d\t%d\t%d\t%s\t%.9e\n"
                    % (Z, N, Zd, Nd, "\t".join(str(c) for c in dec), w))
    return len(rows), stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--astrophomes",
                    default=os.environ.get("ASTROPHOMES_PATH", _DEFAULT_ASTROPHOMES))
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--no-energy-renorm", action="store_true",
                    help="keep the raw redistribution multiplicities (sum_s <x_s> ~ 0.96 "
                         "at the highest eps_r, i.e. a systematic energy sink)")
    args = ap.parse_args()

    models, phenom = load_astrophomes(args.astrophomes)
    em = models.EmpiricalModel()
    os.makedirs(args.out, exist_ok=True)

    n = write_basis(em, os.path.join(args.out, "basis.txt"))
    print("basis.txt:          %d eps_r points" % n)

    n = write_redistribution(em, os.path.join(args.out, "redistribution.txt"),
                             renormalise=not args.no_energy_renorm)
    print("redistribution.txt: %d rows" % n)

    n, stats = write_fragments(models, phenom, em, os.path.join(args.out, "fragments.txt"))
    print("fragments.txt:      %d channels, %d mothers" % (n, len(stats)))
    if 5626 in stats:
        print("  Fe-56: %d channels, raw sum %.4f, <dA> = %.4f  (reference 5.815)"
              % stats[5626])


if __name__ == "__main__":
    main()
