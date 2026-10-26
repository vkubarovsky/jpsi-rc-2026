"""Single-particle efficiency vs generated momentum from vm_to_ntuple.groovy ntuples
(80 float32 per event), for Mariana's GEMC and ours, identical definitions.

    python3 eff_ntuple.py <outdir> "<label>=<glob>" "<label>=<glob>"

Slots k: 0 e' (scattered), 1 p, 2 e+, 3 e- (decay).  Block 19+8k:
found pid charge status quality px py pz (REC::Particle row via MC::GenMatch).
Efficiency = generated particle has a GenMatch partner in REC::Particle with the
right PID, in the expected detector (status: FT 1000-1999, FD 2000-3999).
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NF = 80


def load(pattern):
    fs = sorted(glob.glob(pattern))
    a = np.concatenate([np.fromfile(f, "<f4").reshape(-1, NF) for f in fs])
    return a, len(fs)


def part(a, k):
    px, py, pz = a[:, 1 + 3 * k], a[:, 2 + 3 * k], a[:, 3 + 3 * k]
    p = np.sqrt(px * px + py * py + pz * pz)
    th = np.degrees(np.arctan2(np.hypot(px, py), pz))
    b = 19 + 8 * k
    return p, th, a[:, b] > 0.5, a[:, b + 1].astype(int), np.abs(a[:, b + 3]).astype(int)


SPEC = [  # key, slot, pid, detector status range, name, p bins, theta slices
    ("eFT", 0, 11, (1000, 2000), "scattered e$^-$ in FT", np.arange(0.0, 6.01, 0.25), [(2.0, 3.0), (3.0, 4.0), (4.0, 5.0)]),
    ("em", 3, 11, (2000, 4000), "decay e$^-$ in FD", np.arange(0.0, 9.01, 0.5), [(5, 15), (15, 25), (25, 35)]),
    ("ep", 2, -11, (2000, 4000), "decay e$^+$ in FD", np.arange(0.0, 9.01, 0.5), [(5, 15), (15, 25), (25, 35)]),
    ("p", 1, 2212, (2000, 4000), "proton in FD", np.arange(0.0, 4.01, 0.25), [(5, 20), (20, 40), (40, 60)])]


def effcurve(ok, x, edges):
    n, _ = np.histogram(x, edges)
    k, _ = np.histogram(x[ok], edges)
    e = np.where(n > 0, k / np.maximum(n, 1), np.nan)
    s = np.sqrt(np.maximum(e * (1 - e), 1.0 / np.maximum(n, 1) ** 2) / np.maximum(n, 1))
    return e, s, n


out = sys.argv[1]
os.makedirs(out, exist_ok=True)
S = {}
txt = []
for arg in sys.argv[2:]:
    lab, pat = arg.split("=", 1)
    S[lab], nf = load(pat)
    txt.append(f"{lab}: {nf} files, {len(S[lab])} events")
cols = dict(zip(S, ("C3", "C0")))


def found(a, k, pid, st):
    p, th, f, rpid, rst = part(a, k)
    return p, th, f & (rpid == pid) & (rst >= st[0]) & (rst < st[1]), f


f1, ax = plt.subplots(2, 4, figsize=(18, 7.5), gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.06))
for j, (key, k, pid, st, name, edges, _) in enumerate(SPEC):
    c = 0.5 * (edges[1:] + edges[:-1])
    E = {}
    for lab, a in S.items():
        p, th, ok, anyrec = found(a, k, pid, st)
        e, s, n = effcurve(ok, p, edges)
        E[lab] = (e, s, n)
        g = n > 100
        ax[0, j].errorbar(c[g], e[g], yerr=s[g], fmt="o-", ms=3, lw=1, color=cols[lab], label=lab)
        txt.append(f"  {lab:24s} {key:4s} eff(PID, detector) = {ok.mean():.4f}   any GenMatch = {anyrec.mean():.4f}   N_gen {len(p)}")
    l0, l1 = list(S)
    (e0, s0, n0), (e1, s1, n1) = E[l0], E[l1]
    g = (n0 > 100) & (n1 > 100) & (e0 > 0.02)
    r = e1 / np.where(e0 > 0, e0, np.nan)
    rs = r * np.sqrt((s0 / np.maximum(e0, 1e-9)) ** 2 + (s1 / np.maximum(e1, 1e-9)) ** 2)
    ax[1, j].errorbar(c[g], r[g], yerr=rs[g], fmt="o", ms=3, color="k")
    ax[1, j].axhline(1, color="gray", lw=0.8); ax[1, j].set_ylim(0.7, 1.3)
    ax[1, j].set_ylabel(f"{l1} / {l0}", fontsize=7); ax[1, j].set_xlabel(r"$p_{gen}$ (GeV)")
    ax[0, j].set_title(name, fontsize=10); ax[0, j].set_ylim(0, 1.0); ax[0, j].grid(alpha=0.3)
    ax[0, j].set_xticklabels([]); ax[1, j].set_xlim(ax[0, j].get_xlim())
ax[0, 0].set_ylabel("efficiency = rec / gen"); ax[0, 0].legend(fontsize=8, loc="upper left")
f1.suptitle(r"J/$\psi\to e^+e^-$, F18 inbending, 45 nA background: efficiency vs generated momentum (integrated over $\theta$)", fontsize=11)
f1.savefig(os.path.join(out, "fig_eff_vs_p.png"), dpi=120, bbox_inches="tight")
f1.savefig(os.path.join(out, "fig_eff_vs_p.pdf"), bbox_inches="tight")

f2, ax = plt.subplots(3, 3, figsize=(15, 10), sharey=True)
for j, (key, k, pid, st, name, edges, sl) in enumerate(SPEC[:3]):
    c = 0.5 * (edges[1:] + edges[:-1])
    for i, (t0, t1) in enumerate(sl):
        for lab, a in S.items():
            p, th, ok, _ = found(a, k, pid, st)
            m = (th >= t0) & (th < t1)
            e, s, n = effcurve(ok[m], p[m], edges)
            g = n > 100
            ax[i, j].errorbar(c[g], e[g], yerr=s[g], fmt="o-", ms=3, lw=1, color=cols[lab], label=lab)
            txt.append(f"  {lab:24s} {key:4s} theta {t0:4.1f}-{t1:4.1f}: eff = {ok[m].mean() if m.any() else float('nan'):.4f}  N {m.sum()}")
        ax[i, j].set_title(f"{name}, $\\theta_{{gen}}$ {t0:g}-{t1:g} deg", fontsize=9)
        ax[i, j].set_ylim(0, 1.0); ax[i, j].grid(alpha=0.3)
    ax[2, j].set_xlabel(r"$p_{gen}$ (GeV)")
ax[0, 0].legend(fontsize=8, loc="upper left")
for i in range(3):
    ax[i, 0].set_ylabel("efficiency")
f2.suptitle(r"Efficiency vs $p_{gen}$ in slices of $\theta_{gen}$ (removes the different angular populations of the two generators)", fontsize=11)
f2.savefig(os.path.join(out, "fig_eff_vs_p_theta.png"), dpi=120, bbox_inches="tight")
f2.savefig(os.path.join(out, "fig_eff_vs_p_theta.pdf"), bbox_inches="tight")
open(os.path.join(out, "eff_summary.txt"), "w").write("\n".join(txt) + "\n")
print("\n".join(txt))
