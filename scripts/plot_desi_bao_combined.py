import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

c_km = 299792.458  # km/s

# ============================================================================
# DESI Y1 DATA (arXiv:2404.03002)
# ============================================================================
DESI_DATA = {
    'BGS':  {'z': 0.30, 'DM_rd': 7.93,  'DM_err': 0.15, 'DH_rd': 20.0,  'DH_err': 0.8},
    'LRG1': {'z': 0.51, 'DM_rd': 13.62, 'DM_err': 0.25, 'DH_rd': 22.3,  'DH_err': 0.5},
    'LRG2': {'z': 0.71, 'DM_rd': 17.86, 'DM_err': 0.33, 'DH_rd': 20.4,  'DH_err': 0.5},
    'LRG3': {'z': 0.93, 'DM_rd': 21.71, 'DM_err': 0.28, 'DH_rd': 17.9,  'DH_err': 0.4},
    'ELG':  {'z': 1.32, 'DM_rd': 27.79, 'DM_err': 0.69, 'DH_rd': 13.8,  'DH_err': 0.4},
    'QSO':  {'z': 1.49, 'DM_rd': 30.69, 'DM_err': 0.80, 'DH_rd': 13.3,  'DH_err': 0.5},
    'LyA':  {'z': 2.33, 'DM_rd': 39.71, 'DM_err': 0.94, 'DH_rd': 8.52,  'DH_err': 0.17},
}


def get_bao_from_class(params, label=""):
    """Run CLASS and extract BAO observables at DESI redshifts."""
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()

    rd = cosmo.rs_drag()
    h = cosmo.h()
    z_values = [DESI_DATA[k]['z'] for k in DESI_DATA.keys()]

    DM_rd, DH_rd = [], []
    for z in z_values:
        DA = cosmo.angular_distance(z)       # Mpc
        H_z = cosmo.Hubble(z) * c_km         # km/s/Mpc
        DM = (1 + z) * DA
        DH = c_km / H_z
        DM_rd.append(DM / rd)
        DH_rd.append(DH / rd)

    if label:
        print(f"{label}: h={h:.4f}, Omega_m={cosmo.Omega_m():.4f}, r_d={rd:.2f} Mpc")

    cosmo.struct_cleanup()
    cosmo.empty()
    return np.array(z_values), np.array(DM_rd), np.array(DH_rd)


def create_combined_figure():
    common = {
        'output': 'mPk',
        'P_k_max_1/Mpc': 3.0,
        'omega_b': 0.02242,
        'n_s': 0.9665,
        'A_s': 2.1e-9,
        'tau_reio': 0.054,
    }

    # 1. Planck LCDM
    lcdm_params = {**common, 'h': 0.6736, 'omega_cdm': 0.1200}

    # 2. Holographic (LCDM background with Omega_m ~ 0.318)
    holo_params = {**common, 'h': 0.6784, 'omega_cdm': 0.1193}

    # 3. DESI-preferred w0waCDM
    desi_params = {
        **common,
        'h': 0.68, 'omega_cdm': 0.118,
        'Omega_Lambda': 0.0, 'Omega_k': 0.0,
        'w0_fld': -0.75, 'wa_fld': -0.85,
    }

    z_plot, DM_lcdm, DH_lcdm = get_bao_from_class(lcdm_params, "Planck LCDM")
    _, DM_holo, DH_holo = get_bao_from_class(holo_params, "Holographic")
    _, DM_desi_model, DH_desi_model = get_bao_from_class(desi_params, "DESI w0waCDM")

    DM_data = np.array([DESI_DATA[k]['DM_rd'] for k in DESI_DATA.keys()])
    DM_err  = np.array([DESI_DATA[k]['DM_err'] for k in DESI_DATA.keys()])
    DH_data = np.array([DESI_DATA[k]['DH_rd'] for k in DESI_DATA.keys()])
    DH_err  = np.array([DESI_DATA[k]['DH_err'] for k in DESI_DATA.keys()])

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # --- (a) D_M/r_d ---
    ax = axes[0, 0]
    ax.errorbar(z_plot, DM_data, yerr=DM_err, fmt='ko', ms=8,
                capsize=4, label='DESI Y1 data', zorder=5)
    ax.plot(z_plot, DM_lcdm, 'b-', lw=2.5, label=r'Planck $\Lambda$CDM')
    ax.plot(z_plot, DM_holo, 'r--', lw=2.5, label=r'Holographic ($\Omega_m=0.318$)')
    ax.plot(z_plot, DM_desi_model, 'g:', lw=2.5, label=r'DESI preferred $w_0w_a$CDM')
    ax.set_xlabel('Redshift $z$', fontsize=13)
    ax.set_ylabel(r'$D_M(z)/r_d$', fontsize=13)
    ax.legend(fontsize=10, loc='upper left')
    ax.set_title(r'(a) Comoving Angular Diameter Distance', fontsize=13)
    ax.grid(True, alpha=0.3)

    # --- (b) D_H/r_d ---
    ax = axes[0, 1]
    ax.errorbar(z_plot, DH_data, yerr=DH_err, fmt='ko', ms=8,
                capsize=4, label='DESI Y1 data', zorder=5)
    ax.plot(z_plot, DH_lcdm, 'b-', lw=2.5, label=r'Planck $\Lambda$CDM')
    ax.plot(z_plot, DH_holo, 'r--', lw=2.5, label=r'Holographic ($\Omega_m=0.318$)')
    ax.plot(z_plot, DH_desi_model, 'g:', lw=2.5, label=r'DESI preferred $w_0w_a$CDM')
    ax.set_xlabel('Redshift $z$', fontsize=13)
    ax.set_ylabel(r'$D_H(z)/r_d$', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title(r'(b) Hubble Distance', fontsize=13)
    ax.grid(True, alpha=0.3)

    # --- (c) D_M residuals ---
    ax = axes[1, 0]
    res_lcdm = 100 * (DM_lcdm - DM_data) / DM_data
    res_holo = 100 * (DM_holo - DM_data) / DM_data
    res_desi = 100 * (DM_desi_model - DM_data) / DM_data
    res_err  = 100 * DM_err / DM_data

    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.fill_between(z_plot, -res_err, res_err, alpha=0.25, color='gray', label=r'DESI 1$\sigma$')
    ax.plot(z_plot, res_lcdm, 'bs-', ms=9, lw=2, label=r'Planck $\Lambda$CDM')
    ax.plot(z_plot, res_holo, 'r^--', ms=9, lw=2, label=r'Holographic')
    ax.plot(z_plot, res_desi, 'go:', ms=9, lw=2, label=r'DESI $w_0w_a$CDM')
    ax.set_xlabel('Redshift $z$', fontsize=13)
    ax.set_ylabel(r'$\Delta D_M / D_M^{\rm DESI}$ [%]', fontsize=13)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_title(r'(c) $D_M/r_d$ Residuals', fontsize=13)
    ax.set_ylim(-6, 6)
    ax.grid(True, alpha=0.3)

    # --- (d) D_H residuals ---
    ax = axes[1, 1]
    res_lcdm = 100 * (DH_lcdm - DH_data) / DH_data
    res_holo = 100 * (DH_holo - DH_data) / DH_data
    res_desi = 100 * (DH_desi_model - DH_data) / DH_data
    res_err  = 100 * DH_err / DH_data

    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.fill_between(z_plot, -res_err, res_err, alpha=0.25, color='gray', label=r'DESI 1$\sigma$')
    ax.plot(z_plot, res_lcdm, 'bs-', ms=9, lw=2, label=r'Planck $\Lambda$CDM')
    ax.plot(z_plot, res_holo, 'r^--', ms=9, lw=2, label=r'Holographic')
    ax.plot(z_plot, res_desi, 'go:', ms=9, lw=2, label=r'DESI $w_0w_a$CDM')
    ax.set_xlabel('Redshift $z$', fontsize=13)
    ax.set_ylabel(r'$\Delta D_H / D_H^{\rm DESI}$ [%]', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title(r'(d) $D_H/r_d$ Residuals', fontsize=13)
    ax.set_ylim(-10, 10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('desi_bao_combined.png', dpi=180, bbox_inches='tight')
    plt.savefig('desi_bao_combined.pdf', bbox_inches='tight')
    print("Saved: desi_bao_combined.png/pdf")


if __name__ == '__main__':
    create_combined_figure()
