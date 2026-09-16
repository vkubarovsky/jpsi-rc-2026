"""Regenerate every computed figure in the J/psi radiative-corrections paper.

Writes PDFs into paper/Figures/.  The Feynman-diagram PDFs in that directory
are drawn by hand and are not touched.

Run:  python3 make_figures.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import rc_vector_mesons as rc
from rc_vector_mesons import (ELECTRON_MASS as ME, MUON_MASS as MMU,
                              JPSI_MASS as MJ)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper', 'Figures')
LW, FS = 3.0, 20


def save(name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + '.pdf')
    plt.savefig(path)
    plt.close()
    print('  wrote', os.path.relpath(path, os.path.dirname(OUT)))


def fig_dNdk0(mv=MJ, kmin=0.001, kmax=1.0, N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    x = np.linspace(kmin, kmax, N)
    plt.plot(x, rc.photon_spectrum(x, mv, ME), color='red', lw=LW, label=r'$J/\psi\to e^+e^-$')
    plt.plot(x, rc.photon_spectrum(x, mv, MMU), color='blue', lw=LW, label=r'$J/\psi\to \mu^+\mu^-$')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax.set_yscale('log')
    plt.xlabel(r'$E_\gamma$, GeV', fontsize=30)
    plt.ylabel(r'$\frac{1}{\Gamma_0}\frac{d\Gamma}{dE_\gamma}$', fontsize=30)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle(r'$\frac{1}{\Gamma_0}\frac{d\Gamma}{dE_\gamma}$', fontsize=22)
    save('dNdk0_e_mu')


def fig_soft_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    de1 = np.linspace(0.000001, eps, N)
    plt.plot(de1, rc.soft_photon_correction(de1, mv, ME), color='blue', lw=LW, label=r'$J/\psi\to e^+e^-$')
    plt.plot(de1, rc.soft_photon_correction(de1, mv, MMU), color='red', lw=LW, label=r'$J/\psi\to \mu^+\mu^-$')
    plt.xlabel(r'$\epsilon, GeV$', fontsize=30)
    plt.ylabel(r'$\delta_{RC} $', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle('Soft Photon Radiation Corrections', fontsize=22)
    save('soft_RC')


def fig_hard_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    emax = np.linspace(eps, kmax, N)
    plt.plot(emax, rc.hard_photon_correction(eps, emax, mv, ME), color='blue', lw=LW, label=r'$J/\psi\to e^+e^-$')
    plt.plot(emax, rc.hard_photon_correction(eps, emax, mv, MMU), color='red', lw=LW, label=r'$J/\psi\to \mu^+\mu^-$')
    plt.xlabel(r'$kmax, GeV$', fontsize=30)
    plt.ylabel(r'$\delta_{RC} $', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle('Hard Photon Radiation Corrections', fontsize=22)
    save('hard_RC')


def fig_HR_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    emax = np.linspace(eps, kmax, N)
    plt.plot(emax, rc.delta_vm(emax, mv, ME), color='blue', lw=LW, label=r'$J/\psi\to e^+e^-$')
    plt.plot(emax, rc.delta_vm(emax, mv, MMU), color='red', lw=LW, label=r'$J/\psi\to \mu^+\mu^-$')
    plt.xlabel(r'$kmax, GeV$', fontsize=30)
    plt.ylabel(r'$\delta_{RC} $', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle('Radiation Corrections', fontsize=22)
    save('HR_RC')


def fig_RC_Mv_VM(eps=0.01, kmax=0.1, N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    Mv = np.linspace(0.5, 3.1, N)
    dH = rc.hard_photon_correction(eps, kmax, Mv, ME)
    dR = rc.soft_photon_correction(eps, Mv, ME)
    plt.plot(Mv, dH, color='blue', lw=LW, label='delta_H')
    plt.plot(Mv, dR, color='red', lw=LW, label='delta_R')
    plt.plot(Mv, dH + dR, color='black', lw=LW, label='delta_R+H')
    plt.xlabel(r'$M_V,~ GeV$', fontsize=30)
    plt.ylabel(r'$\delta_{RC} $', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle(r'$Vector~Mesons~RC~~\epsilon=$' + str(eps) + ' kmax=' + str(kmax), fontsize=22)
    save('RC_Mv_VM')


def fig_dNdk0_HR(x0=0.001, eps=0.001, kmax=0.1):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    dR = rc.soft_photon_correction(eps, MJ, ME)
    dH = rc.hard_photon_correction(eps, kmax, MJ, ME)
    dT = rc.delta_vm(kmax, MJ, ME)

    x1 = np.linspace(0, x0, 1000)
    y1 = np.full(1000, np.abs(1 + dR) / x0)
    plt.vlines(x=x0, ymin=0, ymax=y1[0], color='blue',
               label=r'$1+\delta_{RC}(soft)~~$=' + f'{1 + dR:.3f}')
    plt.plot(x1, y1, color='blue', lw=LW, label=r'$~~~~~~\delta_{RC}(soft)~~$=' + f'{dR:.3f}')
    plt.fill_between(x1, y1, color='lightblue')
    plt.yscale('log')

    x2 = np.linspace(eps, kmax, 1000)
    plt.plot(x2, rc.photon_spectrum(x2, MJ, ME), color='red', lw=LW,
             label=r'$~~~~~~\delta_{RC}(hard)=$' + f'{dH:.3f}')
    plt.fill_between(x2, rc.photon_spectrum(x2, MJ, ME), color='coral')

    plt.vlines(x=x0, ymin=0, ymax=y1[0], color='blue',
               label=r'$\delta_{RC}(soft)+\delta_{RC}(hard)$=' + f'{dT:.3f}')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(r'$E_\gamma,~GeV$', fontsize=30)
    plt.ylabel(r'$dN/dE_\gamma$', fontsize=30)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    plt.suptitle(r'Electrons: $dN/dE_{\gamma}  ~~  \epsilon=$' + str(eps) + ' kmax=' + str(kmax), fontsize=22)
    save('dNdk0_HR' + '_' + 'kmax_' + str(kmax))


def fig_VM_BH(N=1000):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    dE = np.linspace(0.01, 0.5, N)
    plt.plot(dE, rc.delta_RC64(dE, MJ), color='blue', lw=LW, label='Formala 64 Bethe-Heitler')
    plt.plot(dE, rc.soft_photon_correction(dE, MJ, ME), color='red', lw=LW, label='Formala 22 Vector mesons')
    plt.xlabel(r'$\epsilon, GeV$', fontsize=30)
    plt.ylabel(r'$\delta_{RC} $', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=FS)
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    save('RC_VM_BH')


if __name__ == '__main__':
    print('regenerating figures ->', OUT)
    fig_dNdk0(mv=MJ, kmin=0.001, kmax=1.0, N=1000)
    fig_soft_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000)
    fig_hard_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000)
    fig_HR_RC(mv=MJ, eps=0.001, kmax=1.0, N=1000)
    fig_RC_Mv_VM(eps=0.01, kmax=0.1, N=1000)
    fig_dNdk0_HR(0.001, 0.001, 0.1)
    fig_VM_BH(N=1000)
    print('done')
