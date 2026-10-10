"""Extract gen/rec M(e+e-) of the J/psi decay leptons from fastMC TruthMatch .dat
files (Mariana's GEMC, generator without radiation: M_gen = 3.0969).

    python3 extract_mee_res.py <dataset dir .../v1/dat> <out.npz>

Per event kept: both decay leptons reconstructed (status >= 1, det FD).
Columns: Q2 xB W t (generated), M_gen, M_rec, p/theta of e+ and e- (gen, rec),
status of e' (FT) and proton.
"""
import glob
import re
import sys

import numpy as np

ME = 0.51099895e-3


def p4(p, th, ph):
    th, ph = np.radians(th), np.radians(ph)
    return np.array([np.sqrt(p * p + ME * ME), p * np.sin(th) * np.cos(ph), p * np.sin(th) * np.sin(ph), p * np.cos(th)])


def mass(a, b):
    q = a + b
    return np.sqrt(max(q[0] ** 2 - q[1] ** 2 - q[2] ** 2 - q[3] ** 2, 0.0))


def main(d, out):
    files = sorted(f for f in glob.glob(d + "/jpsi_tm_p*_*.dat") if re.search(r"_p\d+_(train|val)\.dat$", f))
    rows = []
    nev = 0
    for fn in files:
        with open(fn) as fh:
            lines = [l for l in fh if not l.startswith("#")]
        i = 0
        while i < len(lines):
            h = lines[i].split(); npart = 4
            if i + npart >= len(lines) or len(h) != 7:
                break                                   # truncated last event
            parts = [lines[i + 1 + j].split() for j in range(npart)]
            i += 1 + npart
            nev += 1
            ept, prt, dec = parts[0], parts[1], parts[2:]
            if any(int(x[0]) < 1 or int(x[2]) != 1 for x in dec):
                continue
            g = [p4(float(x[3]), float(x[4]), float(x[5])) for x in dec]
            r = [p4(float(x[7]), float(x[8]), float(x[9])) for x in dec]
            ep = dec[0] if dec[0][1] == "-11" else dec[1]
            em = dec[1] if dec[0][1] == "-11" else dec[0]
            rows.append([float(h[2]), float(h[3]), float(h[4]), float(h[5]), mass(*g), mass(*r),
                         float(ep[3]), float(ep[4]), float(ep[7]), float(ep[8]),
                         float(em[3]), float(em[4]), float(em[7]), float(em[8]),
                         int(ept[0]), int(prt[0])])
    a = np.array(rows)
    names = "Q2 xB W t Mgen Mrec pgen_ep thgen_ep prec_ep threc_ep pgen_em thgen_em prec_em threc_em st_eft st_p".split()
    np.savez_compressed(out, **{n: a[:, k] for k, n in enumerate(names)}, nev_total=nev, nfiles=len(files))
    print(d, "files", len(files), "events", nev, "both leptons rec", len(a))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
