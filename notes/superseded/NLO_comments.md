# Comments on "Fully accounted NLO electromagnetic radiative corrections to the width of J/psi decay"

V. Kubarovsky — 16 September 2026

Equation and table numbers refer to the current draft (`w_jpsi-4`).
"VPK note" refers to *Radiative Corrections to the J/psi photoproduction
reaction*.

---

## 1. Eq. (1) and Eq. (17): replace `E_V` by `m_V`

Both read `1/(2 E_V)`. Since the whole calculation is done in the J/psi rest
frame, `E_V = m_V` identically, and the text already says so. Writing `m_V`
removes a symbol that never takes any other value and makes the normalisation
of Eq. (1) directly comparable with the standard two-body form.

## 2. Eq. (6): use the factorised form, and add the ultrarelativistic limit

Eq. (6) is correct — it is algebraically identical to Eq. (7) of the VPK note
(verified symbolically, difference exactly zero). But the factorised form is
easier to read and to check dimensionally, because the mass dependence appears
only through dimensionless ratios:

    Gamma_B = (4 pi alpha^2 / f_V^2) (m_V/3) (1 + 2 m_l^2/m_V^2) sqrt(1 - 4 m_l^2/m_V^2)

Please also add the ultrarelativistic limit, Eq. (8) of the VPK note, which is
the form actually used to fix f_V from the measured leptonic width:

    Gamma_B = (4 pi alpha^2 / f_V^2) (m_V/3)

Worth stating explicitly: the correction connecting the two is *second* order
in the mass ratio, not first. Writing x = m_l^2/m_V^2,

    (1 + 2x) sqrt(1 - 4x) = 1 - 6 x^2 + O(x^3)

The linear term cancels. This is why the massless limit is so accurate:
4e-15 for J/psi -> e+e-, 8e-6 for J/psi -> mu+mu-. It reaches only 2e-3 for
rho -> mu+mu-.

## 3. Table I: please say what m_V = 1.5 GeV is

Table I and Table II are computed at m_V = 1.5 GeV, which is not the J/psi and
not any physical vector meson. A reader meeting it right after Section II
(which sets m_V = 3.0969 GeV) will be confused.

We understand it to be the last row of Table I of ref. [10], where the tabulated
column is `2E (MeV)` running 500 to 1500 with rho, omega and phi marked — i.e.
2E is the meson mass, since E = m_V/2 at resonance (stated in [10] below their
Eq. (22)). One sentence — "for comparison with ref. [10] we take m_V = 1.5 GeV,
the largest mass tabulated there" — would settle it.

## 4. The comparison with ref. [10] does not hold where it is claimed

The text states that Table I can be "successfully compared" with Table I of
ref. [10] at omega/m_e = 0.01, last line of the first half. At that line the two
disagree, and the disagreement grows with r_cut:

    r_max      ref.[10]     this work      difference
    0.10        -0.087       -0.0866          ~0
    0.25        -0.034       -0.0330          0.001
    0.50        -0.002       +0.000044        0.002
    1.00        +0.021       +0.023891        0.0029

The soft part agrees exactly (both give delta_R = -0.747), so the mass and the
cut-off are being matched correctly. The discrepancy is entirely in the hard
part, and equals 5 alpha / (4 pi) = 0.0029 at r_max = 1.

The reason is that **Eq. (24) of ref. [10] is missing a term**. Their
differential spectrum, Eq. (23), is correct and agrees with Eq. (14) of this
work; the term was lost when that spectrum was integrated analytically. Their
Eq. (26), and their Table I, inherit the omission. Recomputing their table from
their own Eq. (26) reproduces every printed digit, electrons and muons alike,
which confirms the diagnosis.

Suggested replacement for the claim of agreement:

> Our Eq. (14) agrees with expression (23) of [10]. The integrated form, Eq. (24)
> of [10], is however missing a term `(alpha/pi)(3 r_max/2 - r_max^2/4)`, so that
> the totals tabulated in Table I of [10] are low by that amount; at r_max = 1
> the deficit is 5 alpha/(4 pi) = 0.0029. Our results agree with [10] in the
> small-r_max rows, where the missing term is below the quoted precision.

## 5. Please add the closed-form total correction

The paper gives delta_self (8), delta_vert (12), delta_soft (13) and delta_hard
as an unevaluated integral (14), then the exact result (28). There is no compact
closed form for the total ultrarelativistic correction, which is what an
experiment actually applies. We would like it in the paper. The correct
expression — Eq. (24) of ref. [10] with the missing term restored — is

    delta(m_V, omega_cut) = -(4 alpha/pi) [ (1/2 - ln(m_V/m_l)) (ln r + 1/3)
                                          + (1/2)(L_2(r) - pi^2/6)
                                          + 17/72
                                          + R(r)
                                          + r^2/16 - (3/8) r ]

with

    r    = 2 omega_cut / m_V
    R(r) = (r - r^2/4 - 3/4) [ ln( (m_V/m_l) sqrt(1-r) ) - 1/2 ]
    L_2(x) = - int_0^x ln(1-t)/t dt

For the muon channel add the electron-loop vacuum polarization, which is exactly
delta_self of Eq. (8) evaluated with m_e.

The last two terms, `r^2/16 - (3/8) r`, are the ones absent from [10].

Checks this expression passes:

* it reproduces Table I of this work to all six digits at every r_cut, both
  channels, and is independent of omega by construction;
* at r = 1 it gives exactly `delta_self + 3 alpha/(4 pi)`, i.e. delta_self plus
  the FSIC of Eq. (16);
* for J/psi -> e+e- at r = 1 it gives 0.026135, matching the Bardin-Shumeiko
  result of Table III;
* for J/psi -> mu+mu- at r = 1 it gives 0.034017 against the exact 0.034097 —
  the 8e-5 residue is the ultrarelativistic approximation itself, consistent
  with the accuracy discussion in Section IV.

## 6. Reference [9]

    [9] P. Chatagnon, e-Print: 2602.22128 [hep-ex] (2026).

is now published. Please replace with

    [9] P. Chatagnon, Phys. Rev. C 113 (2026) 6, 065203.

---

## Additional point, for consideration

**Eq. (4).** As printed,

    sum_pol eps_alpha(p_V) eps_beta(p_V) = -g_{alpha beta}

is the polarization sum for a *massless* vector. For a massive vector meson it is

    sum_lambda eps_alpha eps_beta^* = -g_{alpha beta} + p_alpha p_beta / m_V^2

The extra term does drop here, because it contracts the lepton current with
p_V = p_1 + p_2 and `ubar(p_1)(p1slash + p2slash)v(p_2) = 0` by the Dirac
equations — so Eq. (5) and everything after it are correct. But the identity as
written is not true on its own, and the surrounding text calls it "the photon"
while the sum is over the three J/psi polarizations. A clause such as "the
p_alpha p_beta term is dropped by current conservation" would fix it.
