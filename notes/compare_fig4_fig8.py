"""NLO Fig.4 and VPK Fig.8 plot the same quantity on different axes.

NLO Fig.4 : delta vs omega_cut/omega_max, omega_max = mV/2 - 2 ml^2/mV
VPK Fig.8 : delta vs k_max in GeV, over 0.001 - 1.0 GeV

This overlays them on the NLO axis and marks how far VPK Fig.8 reaches.
"""
import os, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import rc_vector_mesons as rc
ME, MMU, MJ = rc.ELECTRON_MASS, rc.MUON_MASS, rc.JPSI_MASS

fig, ax = plt.subplots(figsize=(10, 7))
for m, tag, c in ((ME, r'$J/\psi\to e^+e^-$', 'blue'), (MMU, r'$J/\psi\to \mu^+\mu^-$', 'red')):
    wmax = MJ/2 - 2*m**2/MJ
    x = np.linspace(0.002, 1.0, 1000)
    ax.plot(x, rc.delta_vm(x*wmax, MJ, m), color=c, lw=3, label=tag)
    ax.plot(1.0/wmax, rc.delta_vm(1.0, MJ, m), 'o', color=c, ms=9,
            label=tag + r'  $k_{max}=1$ GeV (end of VPK Fig.8)')
ax.axhline(0, color='gray', lw=0.8)
ax.axvspan(1.0/(MJ/2), 1.0, color='gray', alpha=0.12)
ax.text(0.80, -0.26, 'beyond the range\nof VPK Fig. 8', fontsize=12, ha='center', color='dimgray')
ax.set_xlabel(r'$\omega_{cut}/\omega_{max}$', fontsize=20)
ax.set_ylabel(r'$\delta$', fontsize=20)
ax.tick_params(labelsize=14)
ax.grid(True, linestyle='--', linewidth=0.5, color='gray')
ax.legend(fontsize=12, loc='lower right')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'compare_fig4_fig8.pdf')
plt.savefig(out, bbox_inches='tight'); plt.close()
print('wrote', out)
