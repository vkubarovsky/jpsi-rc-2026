"""Vacuum polarisation at the J/psi mass: electron, muon, tau loops (exact one loop,
time-like s = M^2) and a hadronic estimate.  Also the inclusive FSR+vertex
correction 3 alpha/4pi and the check against delta_vm(E) of rc_vector_mesons.py.

    python3 vp_table.py <figdir>

Delta_alpha_l(s) = Re Pi_l(s) = (2 alpha/pi) int_0^1 x(1-x) ln|1 - x(1-x) s/m^2| dx
(renormalised at q^2 = 0; alpha(s) = alpha / (1 - Delta_alpha)); the decay rate
through the photon propagator carries |1 - Pi(s)|^-2 ~ 1 + 2 Delta_alpha.
"""
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALPHA = 1 / 137.035999
ME, MMU, MTAU = 0.51099895e-3, 0.1056583755, 1.77686
MJ = 3.0969


def dalpha_lepton(s, m):
    """Re Pi(s) of one lepton loop, time-like s > 0 (exact, any s/m^2)."""
    X = s / m ** 2
    f = lambda x: x * (1 - x) * math.log(abs(1 - x * (1 - x) * X))  # noqa: E731
    pts = []
    if X > 4:  # threshold crossing: 1 = x(1-x)X at x = (1 -+ beta)/2
        b = math.sqrt(1 - 4 / X)
        pts = [(1 - b) / 2, (1 + b) / 2]
    return 2 * ALPHA / math.pi * quad(f, 0, 1, points=pts or None, limit=400)[0]


def dalpha_lepton_asym(s, m):
    return ALPHA / (3 * math.pi) * (math.log(s / m ** 2) - 5 / 3)


def dalpha_had_spacelike(Q2):
    """Burkhardt-type space-like fit used in exclurad_py/core/virtual.py (Delta_alpha, not 2x)."""
    if Q2 < 1.0:
        a, b, c = -1.345e-9, -2.302e-3, 4.091
    elif Q2 < 64.0:
        a, b, c = -1.512e-3, -2.822e-3, 1.218
    else:
        a, b, c = -1.1344e-3, -3.0680e-3, 0.99992
    return -(a + b * math.log(1 + c * Q2))


s = MJ ** 2
rows = []
for name, m in (("e", ME), ("mu", MMU), ("tau", MTAU)):
    da = dalpha_lepton(s, m)
    rows.append((name, m, s / (4 * m * m), da, dalpha_lepton_asym(s, m) if s > 4 * m * m else float("nan")))
dh = dalpha_had_spacelike(s)
fsr = 3 * ALPHA / (4 * math.pi)

L = math.log(MJ / ME)
# delta_vm(E) of rc_vector_mesons.py (Ehlotzky-Mitter at r = 1)
deltaE = -4 * ALPHA / math.pi * ((1 - 6) / 16 + (0.5 - L) * (1 / 3) + 17 / 72)

out = []
out.append(f"J/psi: M = {MJ} GeV, s = {s:.4f} GeV^2, alpha = 1/137.035999\n")
out.append("loop     m (GeV)     s/4m^2     Delta_alpha (exact)   asymptotic (ln s/m^2 - 5/3)   factor 1+2Delta_alpha")
for name, m, r, da, asy in rows:
    out.append(f"{name:5s} {m:11.6g} {r:11.4g}     {da:.6f}            {asy:.6f}                    {1+2*da:.5f}")
out.append(f"had   (space-like fit at Q2 = s, estimate)   {dh:.6f}                                          {1+2*dh:.5f}")
tot_l = sum(r[3] for r in rows)
out.append(f"\nleptons total Delta_alpha = {tot_l:.6f}  (2x = {2*tot_l:.5f});  with hadrons {tot_l+dh:.6f} (2x = {2*(tot_l+dh):.5f})")
out.append(f"electron / muon loop = {rows[0][3]/rows[1][3]:.3f};  difference e - mu = {rows[0][3]-rows[1][3]:.6f}"
           f"  vs (alpha/3pi) ln(m_mu^2/m_e^2) = {ALPHA/(3*math.pi)*math.log(MMU**2/ME**2):.6f}")
out.append(f"\n3 alpha/4 pi (FSR + vertex, inclusive, KLN) = {fsr:.6f}")
out.append(f"delta_vm(E) (rc_vector_mesons.py)          = {deltaE:.6f}")
out.append(f"3 alpha/4pi + 2 Delta_alpha_e              = {fsr + 2*rows[0][3]:.6f}   -> delta_vm(E) = FSR + electron loop only")
full = (1 + fsr) * (1 + 2 * (tot_l + dh))
full_exact = (1 + fsr) / (1 - (tot_l + dh)) ** 2
out.append(f"\nGamma / Gamma_Born:")
out.append(f"   FSR + vertex only             {1+fsr:.5f}")
out.append(f"   + electron loop (delta_vm)    {(1+fsr)/(1-rows[0][3])**2:.5f}")
out.append(f"   + mu, tau loops               {(1+fsr)/(1-tot_l)**2:.5f}")
out.append(f"   + hadronic (estimate)         {full_exact:.5f}")
out.append(f"B_Born = 5.971 % / factor:  {5.971/(1+fsr):.4f} % (FSR only)   {5.971/full_exact:.4f} % (all, |1-Pi|^-2)")
txt = "\n".join(out)
print(txt)

fd = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(fd, exist_ok=True)
open(os.path.join(fd, "vp_table.txt"), "w").write(txt + "\n")

# figure: Delta_alpha of the three loops vs sqrt(s)
rs = np.linspace(0.3, 5.0, 400)
f, ax = plt.subplots(figsize=(7.5, 4.6))
for (name, m), c in zip((("e", ME), ("mu", MMU), ("tau", MTAU)), ("C0", "C1", "C2")):
    ax.plot(rs, [1e3 * dalpha_lepton(x * x, m) for x in rs], color=c, label=fr"$\{name}$ loop" if name != "e" else "e loop")
ax.plot(rs, [1e3 * dalpha_had_spacelike(x * x) for x in rs], "k--", lw=1, label="hadronic, space-like fit at $Q^2=s$ (estimate)")
ax.axvline(MJ, color="gray", ls=":", lw=1); ax.text(MJ + 0.04, 13, r"$M_{J/\psi}$", color="gray")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$\sqrt{s}$ (GeV)"); ax.set_ylabel(r"Re $\Pi(s)$ = $\Delta\alpha$  ($\times 10^{3}$)")
ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0)); ax.set_title("Vacuum polarisation, time-like, one loop", fontsize=10)
f.savefig(os.path.join(fd, "fig_vp_loops.pdf"), bbox_inches="tight")
f.savefig(os.path.join(fd, "fig_vp_loops.png"), dpi=120, bbox_inches="tight")
