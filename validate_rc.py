"""Independent verification of rc_vector_mesons.py.

The formulas in rc_vector_mesons.py are closed-form integrals.  Here they are
checked against a reference built from scratch by numerical quadrature of the
ultrarelativistic expressions of

    I. Akushevich, A. Ilyichev, V. Kubarovsky, V. Zykunov,
    "Fully accounted NLO electromagnetic radiative corrections to the width
     of J/psi decay",

namely delta = delta_self + delta_vert + delta_soft + delta_hard with

    delta_self  = (2a/pi) (L_f/3 - 5/9)                         their Eq.(8)
    delta_hard  = (a/pi) int du/u [L + ln(1-u) - 1][1 + (1-u)^2] their Eq.(14)

and delta_vert + delta_soft taken with the photon mass cancelled analytically.

Two closed-form anchors are also checked:

  * the total is independent of the soft/hard separator eps
  * over the full phase space  delta(rmax=1) = delta_self + 3a/(4pi),
    the mass- and channel-independent Final-State Inclusive Correction

NOTE on the hard part.  delta_H alone is *not* a physical quantity: the
analytic expression drops terms of order eps at the lower integration limit,
which cancel against terms dropped in delta_S.  So delta_H is compared to the
reference by verifying that the residual vanishes linearly in eps, while the
total is compared strictly.

Run:  python3 validate_rc.py
"""
import sys
import numpy as np
from scipy.integrate import quad

import rc_vector_mesons as rc
from rc_vector_mesons import ALPHA, PI, ELECTRON_MASS as ME, MUON_MASS as MMU, JPSI_MASS as MJ

FAILS = 0


def check(label, got, want, tol):
    global FAILS
    ok = abs(got - want) < tol
    FAILS += (not ok)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:44s} {got:+.9f}  ref {want:+.9f}  d={got - want:+.2e}")


# ------------------------------------------------------------------ reference
def d_self(mV, mf):
    return (2 * ALPHA / PI) * (np.log(mV ** 2 / mf ** 2) / 3 - 5. / 9.)


def d_vert_soft(mV, om, ml):
    L = np.log(mV ** 2 / ml ** 2)
    return (ALPHA / PI) * (1.5 * L - 2 + PI ** 2 / 3) + (2 * ALPHA / PI) * (L - 1) * np.log(2 * om / mV)


def d_hard(mV, om, omcut, ml):
    L = np.log(mV ** 2 / ml ** 2)
    f = lambda u: (1 / u) * (L + np.log(1 - u) - 1) * (1 + (1 - u) ** 2)
    return (ALPHA / PI) * quad(f, 2 * om / mV, 2 * omcut / mV, limit=800)[0]


def ura_total(mV, kmax, ml):
    om = 1e-12 * ml
    d = d_self(mV, ml) + d_vert_soft(mV, om, ml) + d_hard(mV, om, kmax, ml)
    if ml > 0.1:
        d += d_self(mV, ME)
    return d


KMAXS = [0.02, 0.05, 0.1, 0.3, 0.7, 1.2, 1.5]
CHANNELS = ((ME, 'e+e-'), (MMU, 'mu+mu-'))

print("=" * 80)
print(f"1. TOTAL correction delta(kmax), J/psi  mV = {MJ} GeV   [strict, tol 1e-8]")
print("=" * 80)
for ml, tag in CHANNELS:
    print(f"-- {tag}")
    for kmax in KMAXS:
        check(f"delta_vm  kmax={kmax}", rc.delta_vm(kmax, MJ, ml), ura_total(MJ, kmax, ml), 1e-8)

print()
print("=" * 80)
print("2. ANCHOR  full phase space:  delta(rmax=1) == delta_self + 3*alpha/(4*pi)")
print("=" * 80)
FSIC = 0.75 * ALPHA / PI
for ml, tag in CHANNELS:
    want = d_self(MJ, ml) + FSIC + (d_self(MJ, ME) if ml > 0.1 else 0.0)
    check(f"{tag:7s} delta(kmax=mV/2)", rc.delta_vm(MJ / 2, MJ, ml), want, 1e-9)
print(f"       FSIC = 3*alpha/(4*pi) = {FSIC:.6f}  (mass- and channel-independent)")

print()
print("=" * 80)
print("3. ANCHOR  eps-independence of delta_S + delta_H  (kmax = 0.1 GeV)")
print("=" * 80)
for ml, tag in CHANNELS:
    vals = [rc.soft_photon_correction(e, MJ, ml) + rc.hard_photon_correction(e, 0.1, MJ, ml)
            for e in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)]
    spread = max(vals) - min(vals)
    ok = spread <= 1e-9
    FAILS += (not ok)
    print(f"  {'ok  ' if ok else 'FAIL'} {tag:7s} delta = {vals[0]:+.9f}   spread over eps = {spread:.1e}")

print()
print("=" * 80)
print("4. HARD part: residual vs reference must vanish linearly in eps")
print("=" * 80)
for ml, tag in CHANNELS:
    print(f"-- {tag}  (kmax = 0.1 GeV)")
    prev = None
    for eps in (1e-3, 1e-4, 1e-5, 1e-6):
        diff = rc.hard_photon_correction(eps, 0.1, MJ, ml) - d_hard(MJ, eps, 0.1, ml)
        ratio = '' if prev is None else f"   ratio to previous = {prev / diff:6.2f}  (expect 10)"
        print(f"       eps={eps:<8g} residual = {diff:+.3e}{ratio}")
        prev = diff
    # residual must be below the analytic O(eps) bound, with margin
    eps = 1e-5
    r0 = 2 * eps / MJ
    L = np.log(MJ ** 2 / ml ** 2)
    bound = 2.0 * (ALPHA / PI) * ((L - 1) * 2 * r0 + 2 * r0)
    diff = abs(rc.hard_photon_correction(eps, 0.1, MJ, ml) - d_hard(MJ, eps, 0.1, ml))
    ok = diff < bound
    FAILS += (not ok)
    print(f"  {'ok  ' if ok else 'FAIL'} residual {diff:.2e} < O(eps) bound {bound:.2e}")

print()
print("=" * 80)
print("5. SPECTRUM consistency: int dGamma_H/dk over (eps,kmax) == delta_H, up to O(eps)")
print("=" * 80)
for ml, tag in CHANNELS:
    for kmax in (0.1, 0.7, 1.5):
        eps = 1e-6
        num = quad(lambda k: rc.photon_spectrum(k, MJ, ml), eps, kmax, limit=800)[0]
        check(f"{tag:7s} kmax={kmax}", num, rc.hard_photon_correction(eps, kmax, MJ, ml), 1e-6)

print()
print("=" * 80)
print("RESULT:", "ALL CHECKS PASSED" if FAILS == 0 else f"{FAILS} FAILURE(S)")
print("=" * 80)
sys.exit(1 if FAILS else 0)
