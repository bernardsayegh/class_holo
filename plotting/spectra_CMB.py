import sys
sys.path.insert(0, 'python')
from classy import Class
from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt

THETA_S_TARGET = 1.040423
BASE_PARAMS = {
    'output': 'tCl,lCl,mPk', 'lensing': 'yes', 'l_max_scalars': 2600,
    'P_k_max_1/Mpc': 10.0, 'omega_b': 0.02242, 'omega_cdm': 0.11933,
    'n_s': 0.9665, 'A_s': 2.1e-9, 'tau_reio': 0.054, 'N_ur': 2.0328,
    'N_ncdm': 1, 'm_ncdm': 0.06, 'Omega_k': 0.0,
}

def get_cosmo(h, beta):
    cosmo = Class()
    params = dict(BASE_PARAMS)
    params['h'] = h
    if beta > 0:
        params.update({
            'interaction_beta': beta,
            'interaction_ieff_type': 4,
            'f_clust': 0.0,
        })
    cosmo.set(params)
    cosmo.compute()
    return cosmo

def solve_h(beta):
    def f(h):
        cosmo = get_cosmo(h, beta)
        theta = cosmo.get_current_derived_parameters(['100*theta_s'])['100*theta_s']
        cosmo.struct_cleanup(); cosmo.empty()
        return theta - THETA_S_TARGET
    return brentq(f, 0.55, 0.80, xtol=1e-6)

h_lcdm = solve_h(0.0); h_holo = solve_h(1.0/12.0)
print(f"h_lcdm = {h_lcdm:.6f}, h_holo = {h_holo:.6f}")

lcdm = get_cosmo(h_lcdm, 0.0); holo = get_cosmo(h_holo, 1.0/12.0)

sig8_lcdm, Om_lcdm = lcdm.sigma8(), lcdm.Omega_m()
sig8_holo, Om_holo = holo.sigma8(), holo.Omega_m()
S8_lcdm = sig8_lcdm * np.sqrt(Om_lcdm / 0.3)
S8_holo = sig8_holo * np.sqrt(Om_holo / 0.3)
print(f"LCDM:  sigma8={sig8_lcdm:.5f}, Omega_m={Om_lcdm:.5f}, S8={S8_lcdm:.5f}")
print(f"Holo:  sigma8={sig8_holo:.5f}, Omega_m={Om_holo:.5f}, S8={S8_holo:.5f}")

cls_lcdm = lcdm.lensed_cl(2500)
cls_holo = holo.lensed_cl(2500)
ell = cls_lcdm['ell'][2:]
factor = ell * (ell + 1) / (2 * np.pi) * (2.7255e6)**2
tt_lcdm = cls_lcdm['tt'][2:] * factor
tt_holo = cls_holo['tt'][2:] * factor
cmb_ratio = tt_holo / tt_lcdm

fig, (ax2, ax2r) = plt.subplots(2, 1, figsize=(8, 8), height_ratios=[3, 1],
                                  sharex=True, gridspec_kw={'hspace': 0.05})
ax2.plot(ell, tt_lcdm, 'k--', lw=1.5, label=rf'$\Lambda$CDM ($S_8={S8_lcdm:.3f}$)')
ax2.plot(ell, tt_holo, 'b-', lw=1.5, alpha=0.8,
         label=rf'Holographic $\beta=1/12$ ($S_8={S8_holo:.3f}$)')
ax2.set_ylabel(r'$D_\ell^{TT}$ [$\mu$K$^2$]', fontsize=13)
ax2.set_xlim(2, 2500); ax2.set_ylim(0, 6500)
ax2.legend(loc='upper right', fontsize=10)
ax2.set_title(r'CMB Temperature Power Spectrum ($\theta_s$-matched)', fontsize=13)
ax2.grid(True, alpha=0.3)
plt.setp(ax2.get_xticklabels(), visible=False)

ax2r.axhline(1.0, color='k', ls='-', lw=1)
ax2r.axhspan(0.99, 1.01, color='gray', alpha=0.2, label=r'$\pm$1% band')
ax2r.plot(ell, cmb_ratio, 'b-', lw=1.5, label=r'Model A / $\Lambda$CDM')
ax2r.set_xlabel(r'Multipole $\ell$', fontsize=13)
ax2r.set_ylabel('Ratio', fontsize=12)
ax2r.set_ylim(0.97, 1.03)
ax2r.legend(loc='lower left', fontsize=9)
ax2r.grid(True, alpha=0.3)

plt.savefig('spectra_CMB.png', dpi=200, bbox_inches='tight')
plt.savefig('spectra_CMB.pdf', dpi=200, bbox_inches='tight')
print("Saved: spectra_CMB.png/pdf")
lcdm.struct_cleanup(); lcdm.empty()
holo.struct_cleanup(); holo.empty()