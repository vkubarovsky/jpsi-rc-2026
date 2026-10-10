"""M(ee) detector response vs W and t from GEMC (generator without radiation).

    python3 analyze_mee_res.py <outdir> mee_*.npz

dM = M_rec - M_gen per event (both decay leptons reconstructed in the FD).
Core: iterative Gaussian fit within [-2 sigma, +2 sigma] (sigma_core, shift);
tail: fraction with dM < shift - 5 sigma_core; also the 68 % half width and the
median.  Bins in W and in -t (generated values).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
D = {}
for f in sys.argv[2:]:
    name = os.path.basename(f)[4:-4]
    z = np.load(f)
    D[name] = {k: z[k] for k in z.files}
    print(f"{name:12s} events {int(z['nev_total']):8d}  both leptons rec {len(z['Mgen']):7d}")


def core(dm):
    dm = dm[np.abs(dm) < 0.3]
    mu, s = np.median(dm), 0.5 * (np.percentile(dm, 84) - np.percentile(dm, 16))
    for _ in range(4):
        sel = np.abs(dm - mu) < 2 * s
        h, e = np.histogram(dm[sel], 40)
        c = 0.5 * (e[1:] + e[:-1])
        g = lambda x, A, m, w: A * np.exp(-0.5 * ((x - m) / w) ** 2)  # noqa: E731
        try:
            (A, mu, s), cov = curve_fit(g, c, h, p0=[h.max(), mu, s], sigma=np.sqrt(np.maximum(h, 1)))
            s = abs(s)
        except RuntimeError:
            break
    err = np.sqrt(np.diag(cov))
    return mu, s, err[1], err[2]


def stats(dm):
    mu, s, emu, es = core(dm)
    return dict(n=len(dm), shift=mu, sigma=s, eshift=emu, esigma=es,
                tail=np.mean(dm < mu - 5 * s), med=np.median(dm),
                w68=0.5 * (np.percentile(dm, 84) - np.percentile(dm, 16)))


groups = {"F18in": [k for k in D if k.startswith("F18in")], "F18out": [k for k in D if k.startswith("F18out")],
          "S18in": [k for k in D if k.startswith("S18in")], "S18out": [k for k in D if k.startswith("S18out")],
          "S19in": [k for k in D if k.startswith("S19")]}
groups = {g: v for g, v in groups.items() if v}


def cat(keys, col):
    return np.concatenate([D[k][col] for k in keys])


lines = []
lines.append("dataset        N_rec   shift(MeV)  sigma_core(MeV)  w68(MeV)  median(MeV)  tail(<-5sig)")
for g, ks in list(groups.items()) + [("ALL", list(D))]:
    dm = cat(ks, "Mrec") - cat(ks, "Mgen")
    s = stats(dm)
    lines.append(f"{g:12s} {s['n']:7d}   {s['shift']*1e3:7.2f}      {s['sigma']*1e3:6.2f}+-{s['esigma']*1e3:.2f}    {s['w68']*1e3:7.1f}   {s['med']*1e3:7.2f}     {s['tail']:.3f}")

Wb = np.array([4.0, 4.2, 4.3, 4.4, 4.5, 4.6])
tb = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0])
res = {}
for g, ks in groups.items():
    W, t, dm = cat(ks, "W"), -cat(ks, "t"), cat(ks, "Mrec") - cat(ks, "Mgen")
    res[g] = {}
    for var, b, x in (("W", Wb, W), ("t", tb, t)):
        rows = []
        for a, c in zip(b[:-1], b[1:]):
            sel = (x >= a) & (x < c)
            if sel.sum() < 300:
                rows.append(None); continue
            rows.append(stats(dm[sel]))
        res[g][var] = rows
        lines.append(f"\n{g}  vs {var}:" + ("   bin      N     shift  sigma_core  w68   tail" ))
        for (a, c), r in zip(zip(b[:-1], b[1:]), rows):
            if r is None:
                lines.append(f"  {a:.2f}-{c:.2f}   <300 events"); continue
            lines.append(f"  {a:.2f}-{c:.2f} {r['n']:7d}  {r['shift']*1e3:6.2f}  {r['sigma']*1e3:6.2f}+-{r['esigma']*1e3:.2f}  {r['w68']*1e3:6.1f}  {r['tail']:.3f}")
txt = "\n".join(lines)
open(os.path.join(out, "mee_resolution_summary.txt"), "w").write(txt + "\n")
print(txt)

# figure 1: dM distributions, all data, log
f, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for g, ks in groups.items():
    dm = cat(ks, "Mrec") - cat(ks, "Mgen")
    for a, rng in ((ax[0], (-0.15, 0.1)), (ax[1], (-1.5, 0.15))):
        h, e = np.histogram(dm, 250, range=rng)
        a.step(e[:-1], h / h.sum() / np.diff(e), where="post", label=f"{g} ({len(dm)})", lw=1)
ax[0].set_title("core", fontsize=9); ax[1].set_title("full range: low-side tail", fontsize=9)
for a in ax:
    a.set_yscale("log"); a.set_xlabel(r"$M_{rec}(e^+e^-) - M_{gen}(e^+e^-)$ (GeV)"); a.set_ylabel("probability / GeV")
ax[0].legend(fontsize=8)
f.suptitle("J/psi -> e+e-, GEMC (generator without radiation), both decay leptons in the FD", fontsize=10)
f.savefig(os.path.join(out, "fig_dM_by_dataset.png"), dpi=120, bbox_inches="tight")

# figure 2: shift, sigma, tail vs W and t
f, ax = plt.subplots(3, 2, figsize=(12, 10), sharex="col")
for col, (var, b, lab) in enumerate((("W", Wb, "W (GeV)"), ("t", tb, r"$-t$ (GeV$^2$)"))):
    c = 0.5 * (b[1:] + b[:-1])
    for i, (g, R) in enumerate(res.items()):
        rows = R[var]
        ok = np.array([r is not None for r in rows])
        cc = c[ok] + (i - 2) * 0.006 * (b[-1] - b[0])
        get = lambda k: np.array([r[k] for r in rows if r is not None])  # noqa: E731
        ax[0, col].errorbar(cc, get("sigma") * 1e3, yerr=get("esigma") * 1e3, fmt="o-", ms=4, label=g)
        ax[1, col].errorbar(cc, get("shift") * 1e3, yerr=get("eshift") * 1e3, fmt="o-", ms=4)
        ax[2, col].plot(cc, get("tail"), "o-", ms=4)
    ax[0, col].set_ylabel(r"$\sigma_{core}$ (MeV)"); ax[1, col].set_ylabel("core shift (MeV)")
    ax[2, col].set_ylabel(r"tail fraction ($\Delta M < $ shift $- 5\sigma$)"); ax[2, col].set_xlabel(lab)
    for a in ax[:, col]:
        for x in b:
            a.axvline(x, color="0.9", lw=0.6, zorder=0)
ax[0, 0].legend(fontsize=8)
f.suptitle(r"M(ee) detector response vs generated W and $-t$ (GEMC, Mariana's datasets)", fontsize=11)
f.savefig(os.path.join(out, "fig_res_vs_W_t.png"), dpi=120, bbox_inches="tight")
print("figures in", out)
