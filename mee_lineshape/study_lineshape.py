"""How to fit M(ee) of J/psi -> e+e-(gamma): studies behind the report.

    python3 study_lineshape.py <generator LUND dir (uncut Born)> <mee_born.npy (production, with cuts)> <outdir>

Uses dNdMee_v2.cpp through libdNdMee.so (the same code a C++ user would call)
and a faithful Python port of the original dNdMee.cpp for comparison.

Pseudo-data = generator-level M(ee) (MC truth) folded event by event with a
Gaussian of width sigma_true: an idealised detector, so that the fit has a
known answer (M = 3.0969, sigma = sigma_true).  GEMC replaces this later.
"""
import ctypes
import glob
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = ctypes.CDLL(os.path.join(HERE, "libdNdMee.so"))
LIB.dNdMee_full.restype = ctypes.c_double
LIB.dNdMee_full.argtypes = [ctypes.c_double] * 6
LIB.peak_weight.restype = ctypes.c_double
LIB.peak_weight.argtypes = [ctypes.c_double] * 3
ME = 0.51099895e-3
MJ = 3.0969
G4 = np.array([1, -1, -1, -1.0])

# ---------------------------------------------------------------- models
def v2_binned(edges, N, M, sigma, br=-1.0, eps=1e-3, nsub=3):
    lo, hi = edges[:-1], edges[1:]
    out = np.empty(len(lo))
    for i, (a, b) in enumerate(zip(lo, hi)):
        xs = a + (b - a) * (np.arange(nsub) + 0.5) / nsub
        out[i] = np.mean([LIB.dNdMee_full(x, M, sigma, br, eps, ME) for x in xs]) * (b - a)
    return N * out


# faithful port of the ORIGINAL dNdMee.cpp (x = M - mee as photon energy,
# first-order spectrum with fixed Mj, Norm = 2.0707, br free)
_A, _M0, _E0 = 1 / 137.036, 3.096916, 3.096916 / 2


def _dNdk0_orig(k):
    x = k / _E0
    ok = (k >= 1e-5) & (x <= 0.9999)
    out = np.zeros_like(k)
    out[ok] = (2 * _A / np.pi) * (2 * np.log(2 * _E0 / 0.000511 * np.sqrt(1 - x[ok])) - 1) / k[ok] * (1 - x[ok] + x[ok] ** 2 / 2)
    return out


def orig_density(mee, M, sigma, br, eps=1e-3, Norm=2.070655621874216):
    x = M - mee
    g = lambda z: np.exp(-z * z / (2 * sigma * sigma)) / (np.sqrt(2 * np.pi) * sigma)  # noqa: E731
    out = (1 - br) * g(x)
    xmax = np.maximum(eps, x + 3.5 * sigma)
    for i in range(len(mee)):
        if xmax[i] <= eps or x[i] >= 1.5:
            continue
        t = np.linspace(eps, xmax[i], 600)
        out[i] += np.trapezoid(br * Norm * _dNdk0_orig(t) * g(x[i] - t), t)
    return out


def orig_binned(edges, N, M, sigma, br):
    c = 0.5 * (edges[1:] + edges[:-1])
    return N * orig_density(c, M, sigma, br) * np.diff(edges)


# ---------------------------------------------------------------- fitting
def poisson_fit(n, edges, fn, p0, names):
    def nll(p):
        mu = np.maximum(fn(edges, *p), 1e-300)
        return np.sum(mu - n * np.log(mu))
    r = minimize(nll, p0, method="Nelder-Mead",
                 options=dict(xatol=1e-7, fatol=1e-6, maxiter=4000, maxfev=8000))
    p = r.x
    k = len(p)
    h = np.array([max(abs(v) * 1e-4, 1e-6) for v in p])
    H = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            ei, ej = np.eye(k)[i] * h[i], np.eye(k)[j] * h[j]
            H[i, j] = (nll(p + ei + ej) - nll(p + ei - ej) - nll(p - ei + ej) + nll(p - ei - ej)) / (4 * h[i] * h[j])
    try:
        err = np.sqrt(np.abs(np.diag(np.linalg.inv(H))))
    except np.linalg.LinAlgError:
        err = np.full(k, np.nan)
    mu = fn(edges, *p)
    s = mu > 5
    chi2 = float(np.sum((n[s] - mu[s]) ** 2 / mu[s]))
    return dict(zip(names, p)), dict(zip(names, err)), mu, chi2, int(s.sum() - k)


def load_lund(d):
    m = []
    for f in glob.glob(os.path.join(d, "*.lund")):
        L = [l.split() for l in open(f) if l.strip()]
        i = 0
        while i < len(L):
            k = int(L[i][0]); P = L[i + 1:i + 1 + k]; i += 1 + k
            lep = [np.array([float(p[9]), float(p[6]), float(p[7]), float(p[8])]) for p in P
                   if p[3] in ("11", "-11") and p[4] != "0"]
            q = lep[0] + lep[1]
            m.append(np.sqrt(max(q @ (G4 * q), 0)))
    return np.array(m)


# ---------------------------------------------------------------- studies
def main(lunddir, prodnpy, outdir):
    os.makedirs(outdir, exist_ok=True)
    truth = load_lund(lunddir)
    prod = np.load(prodnpy)
    rng = np.random.default_rng(1)
    R = {"n_truth": len(truth), "n_prod": len(prod)}
    edges = np.arange(2.60, 3.30001, 0.005)

    def smear(m, s):
        return m + rng.normal(0, s, len(m))

    # S1: pseudo-data sigma_true = 30 MeV, three fit models
    for sig_true in (0.030, 0.050):
        d = smear(truth, sig_true)
        n, _ = np.histogram(d, edges)
        S = {}
        S["F1_v2_brQED"] = poisson_fit(n, edges, lambda e, N, M, s: v2_binned(e, N, M, s),
                                       [n.sum(), MJ, sig_true * 1.1], ["N", "M", "sigma"])
        S["F2_v2_brfree"] = poisson_fit(n, edges, lambda e, N, M, s, br: v2_binned(e, N, M, s, br),
                                        [n.sum(), MJ, sig_true * 1.1, 0.4], ["N", "M", "sigma", "br"])
        S["F3_orig_brfree"] = poisson_fit(n, edges, lambda e, N, M, s, br: orig_binned(e, N, M, s, br),
                                          [n.sum(), MJ, sig_true * 1.1, 0.4], ["N", "M", "sigma", "br"])
        R[f"smeared_{int(sig_true*1e3)}MeV"] = {k: dict(par=v[0], err=v[1], chi2=v[3], ndf=v[4]) for k, v in S.items()}
        np.save(os.path.join(outdir, f"smeared_{int(sig_true*1e3)}.npy"),
                dict(edges=edges, n=n, mu={k: v[2] for k, v in S.items()}, par={k: v[0] for k, v in S.items()}), allow_pickle=True)
        print(f"sigma_true {sig_true*1e3:.0f} MeV:")
        for k, v in S.items():
            print(f"   {k:16s} " + "  ".join(f"{a}={v[0][a]:.5g}+-{v[1][a]:.2g}" for a in v[0]) + f"  chi2/ndf={v[3]:.0f}/{v[4]}")

    # S2: eps scan on the 30 MeV pseudo-data (F1, br = QED(eps))
    d = smear(truth, 0.030)
    n, _ = np.histogram(d, edges)
    R["eps_scan"] = {}
    for eps in (1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2):
        p, e, mu, chi2, ndf = poisson_fit(n, edges, lambda E, N, M, s: v2_binned(E, N, M, s, -1.0, eps),
                                          [n.sum(), MJ, 0.033], ["N", "M", "sigma"])
        R["eps_scan"][f"{eps*1e3:g}"] = dict(par=p, err=e, chi2=chi2, ndf=ndf,
                                             br_qed=1 - LIB.peak_weight(MJ, eps, ME))
        print(f"   eps {eps*1e3:5.1f} MeV: br_QED {1-LIB.peak_weight(MJ, eps, ME):.4f}  M={p['M']:.5f}  sigma={p['sigma']*1e3:.2f} MeV  chi2/ndf={chi2:.0f}/{ndf}")

    # S3: production Born (lepton cuts) folded with 30 MeV: does acceptance bias F1?
    d = smear(prod, 0.030)
    n, _ = np.histogram(d, edges)
    S = {}
    S["F1_v2_brQED"] = poisson_fit(n, edges, lambda e, N, M, s: v2_binned(e, N, M, s),
                                   [n.sum(), MJ, 0.033], ["N", "M", "sigma"])
    S["F2_v2_brfree"] = poisson_fit(n, edges, lambda e, N, M, s, br: v2_binned(e, N, M, s, br),
                                    [n.sum(), MJ, 0.033, 0.4], ["N", "M", "sigma", "br"])
    R["prod_cuts_30MeV"] = {k: dict(par=v[0], err=v[1], chi2=v[3], ndf=v[4]) for k, v in S.items()}
    np.save(os.path.join(outdir, "prod_30.npy"),
            dict(edges=edges, n=n, mu={k: v[2] for k, v in S.items()}, par={k: v[0] for k, v in S.items()}), allow_pickle=True)
    print("production Born with lepton cuts, sigma_true 30 MeV:")
    for k, v in S.items():
        print(f"   {k:16s} " + "  ".join(f"{a}={v[0][a]:.5g}+-{v[1][a]:.2g}" for a in v[0]) + f"  chi2/ndf={v[3]:.0f}/{v[4]}")

    json.dump(R, open(os.path.join(outdir, "results.json"), "w"), indent=1, default=float)
    np.save(os.path.join(outdir, "truth.npy"), truth)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
