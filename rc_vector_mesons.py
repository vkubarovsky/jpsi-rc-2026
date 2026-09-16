"""First-order QED radiative corrections for V -> l+ l- decay.

Single source of truth for the radiative-correction formulas used in
"Radiative Corrections to the J/psi photoproduction reaction".

Vector-meson corrections follow
    F. Ehlotzky and H. Mitter, Nuovo Cimento 55A (1968) 181,
Bethe-Heitler corrections follow
    Phys. Rev. D 97, 076012 (2018).

Conventions
-----------
mV              vector meson mass; the lepton energy is E = mV/2
particle_mass   ELECTRON_MASS or MUON_MASS, selecting the decay channel
eps             soft/hard photon separator (soft photons have E_gamma < eps)
kmax            largest hard-photon energy retained; rmax = kmax/E = 2*kmax/mV
                and rmax = 1 is the full phase space

Only the sum soft + hard is physical: each piece separately depends on the
unphysical separator eps, and the analytic hard formula is exact only up to
terms of order eps.  `delta_vm` is the eps-independent total.

For the muon channel the electron-loop vacuum polarization is an extra
contribution and is added on top; it is selected automatically by
`particle_mass > 0.1`.
"""
import numpy as np
from scipy.special import spence

# Constants
PI = np.pi
ALPHA = 1. / 137.036
ELECTRON_MASS = 0.000511  # in GeV
MUON_MASS = 0.10565  # in GeV
PROTON_MASS = 0.93827  # in GeV
JPSI_MASS = 3.096916  # M of J/psi


def li2(x):
    """Dilogarithm Li2(x) = -int_0^x ln(1-t)/t dt."""
    return spence(1 - x)


def _log_sqrt_1mr(rmax):
    """log(sqrt(1-rmax)) with the rmax -> 1 endpoint kept finite.

    At rmax = 1 the prefactor (rmax - rmax^2/4 - 3/4) vanishes exactly, so
    clipping the argument turns an 0 * (-inf) = nan into the correct 0.
    """
    return 0.5 * np.log(np.clip(1. - rmax, 1e-300, None))


def _dd(rmax):
    """Non-logarithmic term of the integrated hard-photon correction.

    dd = rmax**2/16 - 3*rmax/8, i.e. the -4*ALPHA/PI-weighted form of
    +(ALPHA/PI) * (3*rmax/2 - rmax**2/4).
    """
    return rmax * (rmax - 6.) / 16.


### F. Ehlotzky, Nuovo Cimento Vol LV A, N 1, (1968)

def photon_spectrum(k, mV, particle_mass):
    """Differential hard-photon spectrum (1/Gamma_0) dGamma/dk.

    k             - photon energy
    mV            - vector meson mass (rho, omega, phi, J/psi)
    particle_mass - electron or muon mass, selecting the decay channel
    """
    E = mV / 2
    x = k / E
    return (2 * ALPHA / PI) * (2 * np.log(2 * E / particle_mass * np.sqrt(1 - x)) - 1) * (1 / k) * (1 - x + x ** 2 / 2)


def soft_photon_correction(eps, mV, particle_mass):
    """Soft photon radiative correction, Eq.22.

    eps - maximum soft photon energy.  Typical choice eps = 0.001 GeV.
    """
    E = mV / 2
    d = -4 * ALPHA / PI * ((np.log(2 * E / particle_mass) - 0.5) * (np.log(E / eps) - 13 / 12) + 17 / 72 - PI * PI / 12)
    if particle_mass > 0.1:
        d = d - 4 * ALPHA / PI * (5 / 18 - 1 / 3 * np.log(2 * E / ELECTRON_MASS))
    return d


def hard_photon_correction(eps, kmax, mV, particle_mass):
    """Hard photon radiative correction, Eq.24.

    eps           - minimum hard photon energy (matches the soft cut-off)
    kmax          - maximum hard photon energy; kmax = mV/2 is full phase space
    mV            - vector meson mass
    particle_mass - electron or muon mass, selecting the decay channel
    """
    E = mV / 2.
    rmax = kmax / E
    R = (rmax - rmax ** 2 / 4. - 3. / 4.) * (np.log(2. * E / particle_mass) + _log_sqrt_1mr(rmax) - 0.5)
    return -4 * ALPHA / PI * (_dd(rmax)
                              + (np.log(2 * E / particle_mass) - 0.5) * (np.log(eps / kmax) + 3. / 4.)
                              + 0.5 * li2(rmax) + R)


def delta_vm(kmax, mV, particle_mass):
    """Total (soft + hard) correction, Eq.26.  Independent of eps.

    kmax          - maximum hard photon energy; kmax = mV/2 is full phase space
    mV            - vector meson mass
    particle_mass - electron or muon mass, selecting the decay channel
    """
    E = mV / 2
    rmax = kmax / E
    R = (rmax - rmax ** 2 / 4 - 3 / 4) * (np.log(2 * E / particle_mass) + _log_sqrt_1mr(rmax) - 0.5)
    d = -4 * ALPHA / PI * (_dd(rmax)
                           + (1 / 2 - np.log(2 * E / particle_mass)) * (np.log(rmax) + 1 / 3)
                           + 1 / 2 * (li2(rmax) - PI ** 2 / 6) + 17 / 72 + R)
    if particle_mass > 0.1:
        d = d - 4 * ALPHA / PI * (5 / 18 - 1 / 3 * np.log(2 * E / ELECTRON_MASS))
    return d


###   Phys. Rev. D 97, 076012 (2018)  -- Bethe-Heitler

def delta_RC64(Delta, Mee):  # Formula 64
    m = ELECTRON_MASS
    beta = np.sqrt(1. - 4 * m * m / (Mee * Mee))
    beta2 = beta * beta
    return -(ALPHA / PI) * ((np.log(4 * Delta * Delta / (m * m)) + np.log((1 - beta) / (1 + beta))) * (
                1 + (1 + beta2) / (2 * beta) * np.log((1 - beta) / (1 + beta))) +
                            (1 - beta) / beta * np.log((1 - beta) / (1 + beta)) + (1 + beta2) / (2 * beta) * (
                                        4 * li2(2 * beta / (1 + beta)) - PI * PI))


def delta_RC64exp(Delta, Mee):
    return np.exp(delta_RC64(Delta, Mee)) - 1


def delta_RC65(Delta, Mee):  # Formula 65
    m = ELECTRON_MASS
    return -(ALPHA / PI) * (np.log(4 * Delta * Delta / (Mee * Mee)) * (1 + np.log(m * m / (Mee * Mee))) - PI * PI / 3)


def delta_RC65exp(Delta, Mee):
    return np.exp(delta_RC65(Delta, Mee)) - 1


def delta_RC66(Delta, Mee):  # Formula 66
    m = ELECTRON_MASS
    beta = np.sqrt(1. - 4 * m * m / (Mee * Mee))
    beta2 = beta * beta
    F = 1. - ALPHA * ALPHA / 3. * np.power(1. + ((1 + beta2) / (2 * beta)) * np.log((1 - beta) / (1 + beta)), 2)
    exp = -ALPHA / PI * (np.log(4 * Delta * Delta / (m * m)) + np.log((1 - beta) / (1 + beta))) * (
                1 + ((1 + beta2) / (2 * beta)) * np.log((1 - beta) / (1 + beta)))
    K = 1 - ALPHA / PI * ((1 - beta) / beta * np.log((1 - beta) / (1 + beta)) + (1 + beta2) / (2 * beta) * (
                4 * li2(2 * beta / (1 + beta)) - PI * PI))
    return F * np.exp(exp) * K - 1.
