# J/psi radiative corrections

First-order QED radiative corrections to the decay `V -> l+ l-`, and the
paper *Radiative Corrections to the J/psi photoproduction reaction*
(V. Kubarovsky).

The corrections to the total width `Gamma(V -> l+l-)` are the radiative
corrections to the cross section of `gamma p -> J/psi p -> l+l- p`, because the
production and the decay factorise.

## Layout

| path | what it is |
|---|---|
| `rc_vector_mesons.py` | the library — every formula lives here, and nowhere else |
| `validate_rc.py` | independent numerical verification of that library |
| `make_figures.py` | regenerates every computed figure in the paper |
| `paper/Jpsi_RC.tex` | paper source |
| `paper/Figures/` | figures (computed ones regenerated, Feynman diagrams hand-drawn) |
| `publish.sh` | copies library + paper + figures into OneDrive |

## Usage

```bash
python3 validate_rc.py     # must print ALL CHECKS PASSED
python3 make_figures.py    # rewrites paper/Figures/*.pdf
cd paper && pdflatex Jpsi_RC.tex && pdflatex Jpsi_RC.tex
./publish.sh               # push results to OneDrive
```

## The library

```python
import rc_vector_mesons as rc

rc.photon_spectrum(k, mV, m)                 # (1/Gamma_0) dGamma_H/dk
rc.soft_photon_correction(eps, mV, m)        # delta_S, photons with E < eps
rc.hard_photon_correction(eps, kmax, mV, m)  # delta_H, eps < E < kmax
rc.delta_vm(kmax, mV, m)                     # delta_S + delta_H, eps-independent
```

`m` is `rc.ELECTRON_MASS` or `rc.MUON_MASS` and selects the decay channel; for
the muon channel the electron-loop vacuum polarization is added automatically.
`kmax = mV/2` (that is `rmax = 1`) is the full phase space.

Only `delta_vm` is physical. `delta_S` and `delta_H` each depend on the
unphysical separator `eps`, and the closed-form `delta_H` is exact only up to
terms of order `eps`; those cancel in the sum.

Bethe-Heitler corrections (`delta_RC64`, `delta_RC65`, `delta_RC66`) follow
Phys. Rev. D **97**, 076012 (2018) and are provided for comparison.

## Validation

`validate_rc.py` checks the closed-form library against a reference built by
numerical quadrature of the ultrarelativistic expressions of Akushevich,
Ilyichev, Kubarovsky and Zykunov, *Fully accounted NLO electromagnetic
radiative corrections to the width of J/psi decay*, plus two analytic anchors:

* the total is independent of `eps` (checked over four decades)
* over the full phase space `delta(rmax=1) = delta_self + 3*alpha/(4*pi)`,
  where `3*alpha/(4*pi)` is the mass- and channel-independent Final-State
  Inclusive Correction

Agreement is at the 1e-14 level for both the `e+e-` and `mu+mu-` channels.

## Reference values

J/psi, `eps = 1 MeV`, electron channel unless noted:

| `kmax` (GeV) | `rmax` | `delta` (e+e-) | `delta` (mu+mu-) |
|---|---|---|---|
| 0.05 | 0.032 | -0.17637 | -0.03400 |
| 0.10 | 0.065 | -0.12605 | -0.01645 |
| 0.30 | 0.194 | -0.05203 | +0.00913 |
| 1.548 | 1.000 | +0.02614 | +0.03402 |

## Sources

* F. Ehlotzky and H. Mitter, Nuovo Cimento **55A** (1968) 181
* Phys. Rev. D **97**, 076012 (2018)
