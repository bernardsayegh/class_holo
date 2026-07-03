#!/usr/bin/env python3
import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

BESTFIT = {
    'output': 'tCl,lCl,mPk', 'lensing': 'yes', 'l_max_scalars': 2600,
    'P_k_max_1/Mpc': 10.0,
    'omega_b': 0.02237, 'omega_cdm': 0.12, 'h': 0.68,
    'A_s': 2.0989032e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
    'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
    'Omega_k': 0.0, 'YHe': 0.2454,
}

# LCDM
params_l = dict(BESTFIT)
lcdm = Class(); lcdm.set(params_l); lcdm.compute()

# Holographic
params_h = dict(BESTFIT)
params_h['interaction_beta'] = 1.0/12.0
params_h['interaction_ieff_type'] = 4
params_h['f_clust'] = 0.0
holo = Class(); holo.set(params_h); holo.compute()

sig8_l, Om_l = lcdm.sigma8(), lcdm.Omega_m()
sig8_h, Om_h = holo.sigma8(), holo.Omega_m()
S8_l = sig8_l * np.sqrt(Om_l / 0.3)
S8_h = sig8_h * np.sqrt(Om_h / 0.3)

print(f"LCDM:  sigma8={sig8_l:.4f}, Omega_m={Om_l:.4f}, S8={S8_l:.3f}")
print(f"Holo:  sigma8={sig8_h:.4f}, Omega_m={Om_h:.4f}, S8={S8_h:.3f}")

h = 0.68
k_vals = np.logspace(-3, np.log10(5.0), 150)
Pk_l = np.array([lcdm.pk_lin(k * h, 0.0) * h**3 for k in k_vals])
Pk_h = np.array([holo.pk_lin(k * h, 0.0) * h**3 for k in k_vals])
ratio = Pk_h / Pk_l

fig, (ax1, ax1r) = plt.subplots(2, 1, figsize=(8, 8), height_ratios=[3, 1],
                                  sharex=True, gridspec_kw={'hspace': 0.05})

ax1.loglog(k_vals, Pk_l, 'k--', lw=2, label=rf'$\Lambda$CDM ($S_8={S8_l:.3f}$)')
ax1.loglog(k_vals, Pk_h, 'b-', lw=2, label=rf'Holographic $\beta=1/12$ ($S_8={S8_h:.3f}$)')
ax1.set_ylabel(r'$P(k)$ [$(h^{-1}$Mpc$)^3$]', fontsize=13)
ax1.set_xlim(1e-3, 5); ax1.set_ylim(10, 3e4)
ax1.legend(loc='lower left', fontsize=10)
ax1.set_title('Matter Power Spectrum (fixed Planck 2018 cosmology)', fontsize=13)
ax1.grid(True, alpha=0.3, which='both')
plt.setp(ax1.get_xticklabels(), visible=False)

ax1r.axhline(1.0, color='k', ls='--', lw=1.5)
ax1r.semilogx(k_vals, ratio, 'b-', lw=2, label=r'Holo/$\Lambda$CDM')
ax1r.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
ax1r.set_ylabel(r'$P_{\rm Holo}/P_{\Lambda{\rm CDM}}$', fontsize=12)
ax1r.set_ylim(0.80, 1.02)
ax1r.legend(loc='lower left', fontsize=9)
ax1r.grid(True, alpha=0.3, which='both')

plt.savefig('spectra_Pk.png', dpi=200, bbox_inches='tight')
plt.savefig('spectra_Pk.pdf', dpi=200, bbox_inches='tight')
print("\nSaved: spectra_Pk.png/pdf")

lcdm.struct_cleanup(); lcdm.empty()
holo.struct_cleanup(); holo.empty()