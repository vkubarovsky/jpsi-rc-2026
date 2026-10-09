"""Figures for the M(ee) lineshape report.  python3 make_figures_mee.py <study dir> <figdir>"""
import ctypes
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from study_lineshape import LIB, ME, MJ, orig_density  # noqa: E402

study, fig = sys.argv[1], sys.argv[2]
os.makedirs(fig, exist_ok=True)
R = json.load(open(os.path.join(study, "results.json")))


def save(f, name):
    f.savefig(os.path.join(fig, name + ".pdf"), bbox_inches="tight")
    f.savefig(os.path.join(fig, name + ".png"), dpi=120, bbox_inches="tight")
    plt.close(f)


# 1. the unfolded spectrum: v2 against the original formula, same br
m = np.linspace(2.0, MJ - 0.0011, 3000)
v2 = np.array([LIB.dNdMee_full(x, MJ, 0.0, -1.0, 1e-3, ME) for x in m])
br = 1 - LIB.peak_weight(MJ, 1e-3, ME)
x = MJ - m
E0 = 3.096916 / 2
xx = x / E0
orig = br * 2.070655621874216 * (2 / 137.036 / np.pi) * (2 * np.log(2 * E0 / 0.000511 * np.sqrt(1 - xx)) - 1) / x * (1 - xx + xx ** 2 / 2)
f, (a1, a2) = plt.subplots(2, 1, figsize=(7.5, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.05))
a1.plot(m, v2, label="dNdMee_v2: exact k, Jacobian, exponentiated (= generator)")
a1.plot(m, orig, "--", label="dNdMee.cpp: x = M - M(ee), first order, Norm 2.0707")
a1.set_yscale("log"); a1.set_ylabel(r"$dN/dM_{ee}$ (1/GeV), photon above 1 MeV")
a1.legend(fontsize=8); a1.set_title(f"unfolded tail, same br = 1 - F(1 MeV) = {br:.4f}", fontsize=9)
a2.plot(m, orig / v2, color="k"); a2.axhline(1, color="gray", lw=0.8)
a2.set_ylabel("original / v2"); a2.set_xlabel(r"$M(e^+e^-)$ (GeV)")
save(f, "fig_spectrum_compare")

# 2. smeared pseudo-data fits: how each fit lies on the MC
for tag, title in (("smeared_30", r"MC truth (200k generator events) smeared event by event with Gauss($\sigma$ = 30 MeV)"),
                   ("smeared_50", r"MC truth (200k generator events) smeared event by event with Gauss($\sigma$ = 50 MeV)"),
                   ("prod_30", r"f18in production Born (1M, lepton cuts) smeared with Gauss($\sigma$ = 30 MeV)")):
    d = np.load(os.path.join(study, tag + ".npy"), allow_pickle=True).item()
    e, n, mu, par = d["edges"], d["n"], d["mu"], d["par"]
    c = 0.5 * (e[1:] + e[:-1]); w = np.diff(e)
    cols = {"F1_v2_brQED": "C0", "F2_v2_brfree": "C2", "F3_orig_brfree": "C3"}
    labs = {"F1_v2_brQED": "F1: v2, br = QED", "F2_v2_brfree": "F2: v2, br free", "F3_orig_brfree": "F3: original dNdMee.cpp, br free"}
    keys = list(mu)
    f = plt.figure(figsize=(12.5, 9.5))
    gs = f.add_gridspec(2 + len(keys), 2, height_ratios=[3.2] + [1.0] * len(keys) + [0.0001], hspace=0.08, wspace=0.62)
    aL = f.add_subplot(gs[0, 0]); aG = f.add_subplot(gs[0, 1])
    for a, logy in ((aL, False), (aG, True)):
        a.errorbar(c, n, yerr=np.sqrt(n), fmt="o", ms=2.5, color="k", lw=0.8, label="pseudo-data (MC)")
        for k in keys:
            a.plot(c, mu[k], color=cols[k], lw=1.2, label=labs[k])
        a.set_xlim(e[0], e[-1]); a.set_xticklabels([])
        if logy:
            a.set_yscale("log"); a.set_ylim(max(n[n > 0].min() * 0.3, 1), n.max() * 3)
            a.set_title("log scale: radiative tail", fontsize=9)
        else:
            a.set_title("linear scale: peak", fontsize=9)
        a.set_ylabel("events / 5 MeV")
    aL.legend(fontsize=8, loc="upper left")
    for r, k in enumerate(keys):
        p = par[k]
        for col in (0, 1):
            a = f.add_subplot(gs[1 + r, col])
            pull = (n - mu[k]) / np.sqrt(np.maximum(mu[k], 1))
            ok = mu[k] > 5
            if col == 0:
                a.errorbar(c[ok], n[ok] / mu[k][ok], yerr=np.sqrt(n[ok]) / mu[k][ok], fmt="o", ms=2, lw=0.6, color=cols[k])
                a.axhline(1, color="gray", lw=0.8); a.set_ylim(0.85, 1.15); a.set_ylabel("MC / fit", fontsize=8)
                chi2 = float(np.sum(pull[ok] ** 2))
                txt = f"{labs[k]}\nM = {p['M']:.5f} GeV, $\\sigma$ = {p['sigma']*1e3:.2f} MeV\nN = {p['N']:.0f}, " + (f"br = {p['br']:.3f}" if "br" in p else "br = 0.3945 (QED)") + f"\n$\\chi^2$/ndf = {chi2:.0f}/{ok.sum()-len(p)}"
                a.text(1.0, 0.5, txt, transform=a.transAxes, fontsize=7.5, va="center", ha="center", color=cols[k],
                       bbox=dict(fc="white", ec=cols[k], lw=0.8))
            else:
                a.bar(c[ok], pull[ok], width=w[ok], color=cols[k], alpha=0.7)
                a.axhline(0, color="gray", lw=0.8); a.set_ylim(-4, 4); a.set_ylabel("pull", fontsize=8)
            a.set_xlim(e[0], e[-1])
            if r < len(keys) - 1:
                a.set_xticklabels([])
            else:
                a.set_xlabel(r"$M(e^+e^-)$ (GeV)")
    f.suptitle(title + "\ntruth: M = 3.0969 GeV, $\\sigma$ = %s MeV, N = %s" % (tag.split("_")[1], "200000" if tag.startswith("smeared") else "1M (all decays)"), fontsize=10)
    save(f, "fig_fit_" + tag)

# 3. eps scan
S = R["eps_scan"]
eps = np.array([float(k) for k in S])
M = np.array([S[k]["par"]["M"] for k in S]); dM = np.array([S[k]["err"]["M"] for k in S])
sg = np.array([S[k]["par"]["sigma"] for k in S]) * 1e3; dsg = np.array([S[k]["err"]["sigma"] for k in S]) * 1e3
brq = np.array([S[k]["br_qed"] for k in S])
f, ax = plt.subplots(1, 3, figsize=(12, 3.6))
ax[0].plot(eps, brq, "o-"); ax[0].set_ylabel(r"br$_{QED}$ = 1 - F($\varepsilon$)")
ax[1].errorbar(eps, (M - MJ) * 1e3, yerr=dM * 1e3, fmt="o-"); ax[1].axhline(0, color="gray", lw=0.8); ax[1].set_ylabel(r"$M_{fit} - M_{J/\psi}$ (MeV)")
ax[2].errorbar(eps, sg, yerr=dsg, fmt="o-"); ax[2].axhline(30, color="gray", lw=0.8); ax[2].set_ylabel(r"$\sigma_{fit}$ (MeV), true 30")
for a in ax:
    a.set_xscale("log"); a.set_xlabel(r"$\varepsilon$ (MeV)")
f.suptitle(r"F1 fit (br fixed to QED at each $\varepsilon$), 30 MeV pseudo-data", fontsize=10)
save(f, "fig_eps_scan")
print("figures in", fig)
