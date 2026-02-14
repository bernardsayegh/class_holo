import sys
sys.path.insert(0, 'python')
from classy import Class
from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt

print("Creating combined P(k) + CMB figure with ratio panels...")

THETA_S_TARGET = 1.040423

BASE_PARAMS = {
    'output': 'tCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': 2600,
    'P_k_max_1/Mpc': 10.0,
    'omega_b': 0.02242,
    'omega_cdm': 0.11933,
    'n_s': 0.9665,
    'A_s': 2.1e-9,
    'tau_reio': 0.054,
    'N_ur': 2.0328,
    'N_ncdm': 1,
    'm_ncdm': 0.06,
    'Omega_k': 0.0,
}

def get_cosmo(h, beta):
    cosmo = Class()
    params = dict(BASE_PARAMS)
    params['h'] = h
    if beta > 0:
        params['interaction_beta'] = beta
        params['f_clust'] = 0.0
        params['interaction_use_particle_horizon'] = 0
        params['super_schwarzschild_correction'] = 'yes'
        params['super_schw_amp'] = 0.0
        params['super_schw_Amap'] = 0.0
        params['super_schw_deltaS'] = 0.03
        params['super_schw_gamma'] = 2.0
    cosmo.set(params)
    cosmo.compute()
    return cosmo

def solve_h(beta):
    def f(h):
        cosmo = get_cosmo(h, beta)
        theta = cosmo.get_current_derived_parameters(['100*theta_s'])['100*theta_s']
        cosmo.struct_cleanup()
        cosmo.empty()
        return theta - THETA_S_TARGET
    return brentq(f, 0.55, 0.80, xtol=1e-6)

def pk_hinvMpc3(cosmo, k_hmpc):
    h = cosmo.h()
    return cosmo.pk_lin(k_hmpc * h, 0.0) * h**3

# Solve for h
print("Solving h for LCDM...")
h_lcdm = solve_h(0.0)
print(f"  h = {h_lcdm:.5f}")

print("Solving h for Model A (beta=1/12)...")
h_holo = solve_h(1.0/12.0)
print(f"  h = {h_holo:.5f}")

# Get cosmologies
lcdm = get_cosmo(h_lcdm, 0.0)
holo = get_cosmo(h_holo, 1.0/12.0)

# Derived params
sig8_lcdm, Om_lcdm = lcdm.sigma8(), lcdm.Omega_m()
sig8_holo, Om_holo = holo.sigma8(), holo.Omega_m()
S8_lcdm = sig8_lcdm * np.sqrt(Om_lcdm / 0.3)
S8_holo = sig8_holo * np.sqrt(Om_holo / 0.3)

print(f"LCDM: sigma8={sig8_lcdm:.4f}, S8={S8_lcdm:.4f}")
print(f"Holo: sigma8={sig8_holo:.4f}, S8={S8_holo:.4f}")

# P(k) data
k_vals = np.logspace(-3, np.log10(5.0), 150)
Pk_lcdm = np.array([pk_hinvMpc3(lcdm, k) for k in k_vals])
Pk_holo = np.array([pk_hinvMpc3(holo, k) for k in k_vals])
pk_ratio = Pk_holo / Pk_lcdm

# CMB data
cls_lcdm = lcdm.lensed_cl(2500)
cls_holo = holo.lensed_cl(2500)
ell = cls_lcdm['ell'][2:]
factor = ell * (ell + 1) / (2 * np.pi) * (2.7255e6)**2
tt_lcdm = cls_lcdm['tt'][2:] * factor
tt_holo = cls_holo['tt'][2:] * factor
cmb_ratio = tt_holo / tt_lcdm

# DES band for P(k) ratio
S8_DES, S8_DES_err = 0.776, 0.017
sig8_des_lo = (S8_DES - S8_DES_err) / np.sqrt(Om_holo / 0.3)
sig8_des_hi = (S8_DES + S8_DES_err) / np.sqrt(Om_holo / 0.3)
ratio_des_lo = (sig8_des_lo / sig8_lcdm) ** 2
ratio_des_hi = (sig8_des_hi / sig8_lcdm) ** 2

# Create figure with GridSpec
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.05, wspace=0.25)

ax1 = fig.add_subplot(gs[0, 0])
ax1r = fig.add_subplot(gs[1, 0], sharex=ax1)
ax2 = fig.add_subplot(gs[0, 1])
ax2r = fig.add_subplot(gs[1, 1], sharex=ax2)

# Panel (a) top: P(k)
ax1.loglog(k_vals, Pk_lcdm, 'k--', lw=2, label=rf'$\Lambda$CDM ($S_8$={S8_lcdm:.3f})')
ax1.loglog(k_vals, Pk_holo, 'b-', lw=2, label=rf'Holographic $\beta=1/12$ ($S_8$={S8_holo:.3f})')
ax1.set_ylabel(r'$P(k)$ [$(h^{-1}$Mpc$)^3$]', fontsize=13)
ax1.set_xlim(1e-3, 5)
ax1.set_ylim(10, 3e4)
ax1.legend(loc='lower left', fontsize=10)
ax1.set_title('(a) Matter Power Spectrum', fontsize=13)
ax1.grid(True, alpha=0.3, which='both')
plt.setp(ax1.get_xticklabels(), visible=False)

# Panel (a) bottom: P(k) ratio
ax1r.axhline(1.0, color='k', ls='--', lw=1.5)
ax1r.axhspan(ratio_des_lo, ratio_des_hi, color='lightgreen', alpha=0.5, zorder=0, label=rf'DES-Y3 $S_8$={S8_DES}$\pm${S8_DES_err}')
ax1r.semilogx(k_vals, pk_ratio, 'b-', lw=2, label=r'Holo/$\Lambda$CDM')
ax1r.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
ax1r.set_ylabel(r'$P_{\rm Holo}/P_{\Lambda{\rm CDM}}$', fontsize=12)
ax1r.set_ylim(0.80, 1.02)
ax1r.legend(loc='lower left', fontsize=9)
ax1r.grid(True, alpha=0.3, which='both')

# Panel (b) top: CMB
ax2.plot(ell, tt_lcdm, 'k--', lw=1.5, label=rf'$\Lambda$CDM ($S_8$={S8_lcdm:.3f})')
ax2.plot(ell, tt_holo, 'b-', lw=1.5, alpha=0.8, label=rf'Holographic $\beta=1/12$ ($S_8$={S8_holo:.3f})')
ax2.set_ylabel(r'$D_\ell^{TT}$ [$\mu$K$^2$]', fontsize=13)
ax2.set_xlim(2, 2500)
ax2.set_ylim(0, 6500)
ax2.legend(loc='upper right', fontsize=10)
ax2.set_title(r'(b) CMB Temperature Power Spectrum', fontsize=13)
ax2.grid(True, alpha=0.3)
plt.setp(ax2.get_xticklabels(), visible=False)

# Panel (b) bottom: CMB ratio
ax2r.axhline(1.0, color='k', ls='-', lw=1)
ax2r.axhspan(0.99, 1.01, color='gray', alpha=0.2, label=r'$\pm$1% band')
ax2r.plot(ell, cmb_ratio, 'b-', lw=1.5, label=r'Model A / $\Lambda$CDM')
ax2r.set_xlabel(r'Multipole $\ell$', fontsize=13)
ax2r.set_ylabel('Ratio', fontsize=12)
ax2r.set_ylim(0.97, 1.03)
ax2r.legend(loc='lower left', fontsize=9)
ax2r.grid(True, alpha=0.3)

plt.savefig('spectra_combined.png', dpi=200, bbox_inches='tight')
plt.savefig('spectra_combined.pdf', dpi=200, bbox_inches='tight')
print("Saved: spectra_combined.png/pdf")

lcdm.struct_cleanup(); lcdm.empty()
holo.struct_cleanup(); holo.empty()
