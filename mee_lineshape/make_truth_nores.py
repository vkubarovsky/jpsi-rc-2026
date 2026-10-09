"""MC truth M(ee) without resolution against the lineshape, in dM = M - M(ee).

    python3 make_truth_nores.py <study dir with truth.npy> <figdir>

Log bins in dM make the eps = 1 MeV threshold of the generator visible: decays
with a photon below eps sit exactly at M (dM = 0), and nothing lies between
dM = 0 and dM(k = eps) = 1.0002 MeV.  v2 has NO free parameter here (N = number
of generated decays, br = QED); the original dNdMee.cpp is shown with its own
best fit (br = 0.3802, N = 195353 from the sigma = 0.1 MeV fit).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from study_lineshape import LIB, ME, MJ, _dNdk0_orig  # noqa: E402

study, figdir = sys.argv[1], sys.argv[2]
truth = np.load(os.path.join(study, "truth.npy"))
NTOT = len(truth)
EPS = 1e-3
dM = MJ - truth
peak = np.abs(dM) < 1e-5
n_peak = int(peak.sum())
Fe = LIB.peak_weight(MJ, EPS, ME)
dm_eps = MJ - np.sqrt(MJ * MJ - 2 * MJ * EPS)
ORIG_N, ORIG_BR, ORIG_NORM = 195353.0, 0.3802, 2.070655621874216


def v2_bin(a, b, nsub=40):
    # expected events of v2 (br = QED, no free parameter) in dM in [a, b]
    u = np.exp(np.linspace(np.log(a), np.log(b), nsub + 1))
    c = np.sqrt(u[1:] * u[:-1])
    f = np.array([LIB.dNdMee_full(MJ - x, MJ, 0.0, -1.0, EPS, ME) for x in c])
    return NTOT * np.sum(f * np.diff(u))


def orig_bin(a, b, nsub=40):
    u = np.exp(np.linspace(np.log(a), np.log(b), nsub + 1))
    c = np.sqrt(u[1:] * u[:-1])
    f = ORIG_BR * ORIG_NORM * _dNdk0_orig(c) * (c >= EPS)
    return ORIG_N * np.sum(f * np.diff(u))


edges = np.logspace(-4, np.log10(1.2), 81)
n, _ = np.histogram(dM[~peak], edges)
mu_v2 = np.array([v2_bin(a, b) for a, b in zip(edges[:-1], edges[1:])])
mu_or = np.array([orig_bin(a, b) for a, b in zip(edges[:-1], edges[1:])])
dl = np.diff(np.log10(edges))
c = np.sqrt(edges[1:] * edges[:-1])

f = plt.figure(figsize=(13, 8.5))
gs = f.add_gridspec(2, 2, height_ratios=[3, 1.3], hspace=0.06, wspace=0.22)
a1 = f.add_subplot(gs[0, 0]); a2 = f.add_subplot(gs[1, 0], sharex=a1)
ok = n > 0
a1.errorbar(c[ok], n[ok] / dl[ok], yerr=np.sqrt(n[ok]) / dl[ok], fmt="o", ms=3, color="k",
            label=f"generator, MC truth ({NTOT} decays, no cuts)")
a1.step(edges, np.append(mu_v2 / dl, np.nan), where="post", color="C0", lw=1.4,
        label="dNdMee_v2: br = QED, N = generated (no free parameter)")
a1.step(edges, np.append(mu_or / dl, np.nan), where="post", color="C3", lw=1.2, ls="--",
        label=f"original dNdMee.cpp, best fit: br = {ORIG_BR}, N = {ORIG_N:.0f}")
a1.axvline(dm_eps, color="gray", ls=":", lw=1)
a1.text(dm_eps * 1.08, 30, r"$k = \varepsilon$ = 1 MeV" "\n" r"($\Delta M$ = 1.0002 MeV)", fontsize=8, color="gray")
a1.text(1.3e-4, 2.0e3,
        f"exactly at M ($\\Delta M$ = 0, photon below 1 MeV):\n"
        f"  generator {n_peak} = {n_peak/NTOT:.4f}\n  QED F(1 MeV) = {Fe:.4f}  ->  {Fe*NTOT:.0f}\n"
        f"between 0 and 1 MeV: {int(((dM > 1e-5) & (dM < dm_eps*0.999)).sum())} events",
        fontsize=8.5, va="bottom", bbox=dict(fc="white", ec="gray", lw=0.6))
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_ylim(10, 3e5)
a1.set_ylabel(r"events / unit $\log_{10}\Delta M$")
a1.legend(fontsize=8, loc="upper right")
a1.set_title(r"tail in $\Delta M = M_{J/\psi} - M(e^+e^-)$, log scale", fontsize=9)
plt.setp(a1.get_xticklabels(), visible=False)
s = mu_v2 > 0
a2.errorbar(c[ok & s], n[ok & s] / mu_v2[ok & s], yerr=np.sqrt(n[ok & s]) / mu_v2[ok & s], fmt="o", ms=3, color="C0", label="MC / v2")
so = mu_or > 0
a2.errorbar(c[ok & so] * 1.03, n[ok & so] / mu_or[ok & so], yerr=np.sqrt(n[ok & so]) / mu_or[ok & so], fmt="s", ms=2.5, color="C3", label="MC / original")
a2.axhline(1, color="gray", lw=0.8); a2.axvline(dm_eps, color="gray", ls=":", lw=1)
a2.set_ylim(0.75, 1.35); a2.set_ylabel("MC / function"); a2.legend(fontsize=8, ncol=2, loc="upper left")
a2.set_xlabel(r"$\Delta M = M - M(e^+e^-)$ (GeV)")
chi_v2 = float(np.sum((n[s] - mu_v2[s]) ** 2 / mu_v2[s]))
a2.text(0.99, 0.08, f"v2: $\\chi^2$/ndf = {chi_v2:.0f}/{s.sum()} (0 free parameters)", transform=a2.transAxes, ha="right", fontsize=8, color="C0")

# right: M(ee) next to the peak, 0.1 MeV bins, linear mass axis
b1 = f.add_subplot(gs[0, 1]); b2 = f.add_subplot(gs[1, 1], sharex=b1)
e2 = np.arange(3.08685, 3.09700, 0.0001)   # M = 3.0969 in the middle of the last bin
n2, _ = np.histogram(truth, e2)
c2 = 0.5 * (e2[1:] + e2[:-1])
mu2 = np.array([NTOT * np.mean([LIB.dNdMee_full(x, MJ, 0.0, -1.0, EPS, ME) for x in np.linspace(a, b, 6)[1:-1]]) * (b - a)
                for a, b in zip(e2[:-1], e2[1:])])
ip = np.searchsorted(e2, MJ) - 1
mu2[ip] += Fe * NTOT
z = n2 > 0
b1.errorbar(c2[z], n2[z], yerr=np.sqrt(n2[z]), fmt="o", ms=2.5, color="k", label="generator")
b1.plot(c2[~z], np.full((~z).sum(), 0.7), "v", ms=4, color="k", mfc="none", label="empty bin")
b1.step(e2, np.append(mu2, np.nan), where="post", color="C0", lw=1.2, label="dNdMee_v2 (delta peak in the last bin)")
b1.axvline(MJ - dm_eps, color="gray", ls=":", lw=1)
b1.text(MJ - dm_eps - 0.0001, 3e4, "M - 1 MeV", ha="right", fontsize=8, color="gray")
b1.set_yscale("log"); b1.set_ylim(0.5, 3e5)
b1.set_ylabel("events / 0.1 MeV"); b1.legend(fontsize=8, loc="upper left")
b1.set_title(r"last 10 MeV below $M_{J/\psi}$, 0.1 MeV bins: the empty 1 MeV gap", fontsize=9)
plt.setp(b1.get_xticklabels(), visible=False)
s2 = mu2 > 5
b2.errorbar(c2[s2], n2[s2] / mu2[s2], yerr=np.sqrt(n2[s2]) / mu2[s2], fmt="o", ms=2.5, color="C0")
b2.axhline(1, color="gray", lw=0.8); b2.set_ylim(0.75, 1.25); b2.set_ylabel("MC / v2")
b2.set_xlabel(r"$M(e^+e^-)$ (GeV)")
b2.ticklabel_format(axis="x", useOffset=False)
f.suptitle(r"J/$\psi\to e^+e^-(\gamma)$, generator level, no resolution: MC truth vs the lineshape", fontsize=11)
for ext, kw in (("pdf", {}), ("png", dict(dpi=120))):
    f.savefig(os.path.join(figdir, "fig_truth_nores." + ext), bbox_inches="tight", **kw)
print("peak", n_peak, n_peak / NTOT, "QED", Fe, "gap", int(((dM > 1e-5) & (dM < dm_eps * 0.999)).sum()), "chi2 v2", chi_v2, s.sum())
