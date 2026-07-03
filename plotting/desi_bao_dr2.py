#!/usr/bin/env python3
"""
DESI DR2 BAO Compatibility Figure
==================================
Data from Table IV of arXiv:2503.14738 (DESI Collaboration 2025).
Run from ~/class_holo_test:  python3 desi_bao_dr2.py
"""
import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

c_km = 299792.458

BGS = {'z': 0.295, 'DV_rd': 7.94, 'DV_err': 0.075}

DESI_DR2 = {
    'LRG1':       {'z': 0.510, 'DM_rd': 13.587, 'DM_err': 0.167, 'DH_rd': 20.980, 'DH_err': 0.425},
    'LRG2':       {'z': 0.706, 'DM_rd': 17.347, 'DM_err': 0.180, 'DH_rd': 19.458, 'DH_err': 0.332},
    'LRG3+ELG1': {'z': 0.934, 'DM_rd': 21.574, 'DM_err': 0.153, 'DH_rd': 17.641, 'DH_err': 0.193},
    'ELG2':       {'z': 1.321, 'DM_rd': 27.605, 'DM_err': 0.320, 'DH_rd': 14.178, 'DH_err': 0.217},
    'QSO':        {'z': 1.484, 'DM_rd': 30.519, 'DM_err': 0.758, 'DH_rd': 12.816, 'DH_err': 0.513},
    'LyA':        {'z': 2.330, 'DM_rd': 38.988, 'DM_err': 0.531, 'DH_rd':  8.632, 'DH_err': 0.101},
}


def get_bao_from_class(params, z_values, label=""):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    rd = cosmo.rs_drag()
    h  = cosmo.h()
    DM_rd, DH_rd, DV_rd = [], [], []
    for z in z_values:
        DA  = cosmo.angular_distance(z)
        H_z = cosmo.Hubble(z) * c_km
        DM  = (1 + z) * DA
        DH  = c_km / H_z
        DV  = (z * DM**2 * DH)**(1./3.)
        DM_rd.append(DM / rd)
        DH_rd.append(DH / rd)
        DV_rd.append(DV / rd)
    if label:
        print(f"{label}: h={h:.4f}, Omega_m={cosmo.Omega_m():.4f}, r_d={rd:.2f} Mpc")
    cosmo.struct_cleanup()
    cosmo.empty()
    return np.array(DM_rd), np.array(DH_rd), np.array(DV_rd)


def create_combined_figure():
    common = {
        'output': 'mPk',
        'P_k_max_1/Mpc': 3.0,
        'omega_b': 0.02237,
        'n_s': 0.9665,
        'A_s': 2.1e-9,
        'tau_reio': 0.054,
        'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,   # Planck-nu convention (v5)
    }
    lcdm_params = {**common, 'h': 0.6736, 'omega_cdm': 0.1200}
    holo_params = {**common, 'h': 0.68, 'omega_cdm': 0.12,
                   'interaction_beta': 1.0/12.0,        # interaction ON: this IS the holographic model
                   'interaction_ieff_type': 4,
                   'f_clust': 0.0}

    # DESI DR2 + CMB + DESY5 headline (Table X, 2503.14738)
    # Self-consistent best-fit: w0=-0.752, wa=-0.86, h=0.683, omega_cdm=0.1197
    desi_params = {
        **common,
        'h': 0.683, 'omega_cdm': 0.1197,
        'omega_b': 0.02243,
        'Omega_Lambda': 0.0, 'Omega_k': 0.0,
        'w0_fld': -0.752, 'wa_fld': -0.86,
    }

    z_bgs = [BGS['z']]
    z_aniso = [DESI_DR2[k]['z'] for k in DESI_DR2.keys()]
    z_all = z_bgs + z_aniso

    DM_lcdm, DH_lcdm, DV_lcdm = get_bao_from_class(lcdm_params, z_all, "Planck LCDM")
    DM_holo, DH_holo, DV_holo = get_bao_from_class(holo_params, z_all, "Holographic")
    DM_desi, DH_desi, DV_desi = get_bao_from_class(desi_params, z_all, "DESI w0waCDM")

    z_aniso = np.array(z_aniso)
    DM_data = np.array([DESI_DR2[k]['DM_rd'] for k in DESI_DR2.keys()])
    DM_err  = np.array([DESI_DR2[k]['DM_err'] for k in DESI_DR2.keys()])
    DH_data = np.array([DESI_DR2[k]['DH_rd'] for k in DESI_DR2.keys()])
    DH_err  = np.array([DESI_DR2[k]['DH_err'] for k in DESI_DR2.keys()])
    DV_bgs_data = BGS['DV_rd']
    DV_bgs_err  = BGS['DV_err']

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    ax = axes[0, 0]
    ax.errorbar(z_aniso, DM_data, yerr=DM_err, fmt='ko', ms=8, capsize=4, label='DESI DR2 data', zorder=5)
    ax.plot(z_aniso, DM_lcdm[1:], 'b-',  lw=2.5, label='Planck LCDM')
    ax.plot(z_aniso, DM_holo[1:], 'r--', lw=2.5, label='Holographic (Om ~ 1/3)')
    ax.plot(z_aniso, DM_desi[1:], 'g:',  lw=2.5, label='w0waCDM (DESI+CMB+DESY5)')
    ax.set_xlabel('Redshift z', fontsize=13)
    ax.set_ylabel('D_M(z)/r_d', fontsize=13)
    ax.legend(fontsize=10, loc='upper left')
    ax.set_title('(a) Comoving Angular Diameter Distance', fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.errorbar(z_aniso, DH_data, yerr=DH_err, fmt='ko', ms=8, capsize=4, label='DESI DR2 data', zorder=5)
    ax.plot(z_aniso, DH_lcdm[1:], 'b-',  lw=2.5, label='Planck LCDM')
    ax.plot(z_aniso, DH_holo[1:], 'r--', lw=2.5, label='Holographic (Om ~ 1/3)')
    ax.plot(z_aniso, DH_desi[1:], 'g:',  lw=2.5, label='w0waCDM (DESI+CMB+DESY5)')
    ax.set_xlabel('Redshift z', fontsize=13)
    ax.set_ylabel('D_H(z)/r_d', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title('(b) Hubble Distance', fontsize=13)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    res_lcdm = 100 * (DM_lcdm[1:] - DM_data) / DM_data
    res_holo = 100 * (DM_holo[1:] - DM_data) / DM_data
    res_desi = 100 * (DM_desi[1:] - DM_data) / DM_data
    res_err  = 100 * DM_err / DM_data
    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.fill_between(z_aniso, -res_err, res_err, alpha=0.25, color='gray', label='DESI DR2 1sig')
    ax.plot(z_aniso, res_lcdm, 'bs-',  ms=9, lw=2, label='Planck LCDM')
    ax.plot(z_aniso, res_holo, 'r^--', ms=9, lw=2, label='Holographic')
    ax.plot(z_aniso, res_desi, 'go:',  ms=9, lw=2, label='DESI w0waCDM')
    ax.set_xlabel('Redshift z', fontsize=13)
    ax.set_ylabel('delta D_M / D_M(DESI) [%]', fontsize=13)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_title('(c) D_M/r_d Residuals', fontsize=13)
    ax.set_ylim(-6, 6)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    res_lcdm = 100 * (DH_lcdm[1:] - DH_data) / DH_data
    res_holo = 100 * (DH_holo[1:] - DH_data) / DH_data
    res_desi = 100 * (DH_desi[1:] - DH_data) / DH_data
    res_err  = 100 * DH_err / DH_data
    ax.axhline(0, color='k', ls='-', lw=0.5)
    ax.fill_between(z_aniso, -res_err, res_err, alpha=0.25, color='gray', label='DESI DR2 1sig')
    ax.plot(z_aniso, res_lcdm, 'bs-',  ms=9, lw=2, label='Planck LCDM')
    ax.plot(z_aniso, res_holo, 'r^--', ms=9, lw=2, label='Holographic')
    ax.plot(z_aniso, res_desi, 'go:',  ms=9, lw=2, label='DESI w0waCDM')
    ax.set_xlabel('Redshift z', fontsize=13)
    ax.set_ylabel('delta D_H / D_H(DESI) [%]', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title('(d) D_H/r_d Residuals', fontsize=13)
    ax.set_ylim(-10, 10)
    ax.grid(True, alpha=0.3)

    plt.suptitle('DESI DR2 BAO Compatibility (arXiv:2503.14738)', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig('desi_bao_dr2.png', dpi=180, bbox_inches='tight')
    plt.savefig('desi_bao_dr2.pdf', bbox_inches='tight')
    print("\nSaved: desi_bao_dr2.png/pdf")

    print("\n--- Chi-squared per tracer (DM + DH) ---")
    models = {
        'Planck LCDM':  (DM_lcdm[1:], DH_lcdm[1:]),
        'Holographic':  (DM_holo[1:], DH_holo[1:]),
        'DESI w0waCDM': (DM_desi[1:], DH_desi[1:]),
    }
    for mname, (dm_mod, dh_mod) in models.items():
        chi2_DM = np.sum(((dm_mod - DM_data) / DM_err)**2)
        chi2_DH = np.sum(((dh_mod - DH_data) / DH_err)**2)
        print(f"  {mname:20s}: chi2_DM={chi2_DM:.2f}, chi2_DH={chi2_DH:.2f}, total={chi2_DM+chi2_DH:.2f}  (12 dof)")

    print(f"\n--- BGS DV/rd at z=0.295 ---")
    print(f"  Data:         {DV_bgs_data:.3f} +/- {DV_bgs_err:.3f}")
    print(f"  Planck LCDM:  {DV_lcdm[0]:.3f}  (pull = {(DV_lcdm[0]-DV_bgs_data)/DV_bgs_err:+.2f} sigma)")
    print(f"  Holographic:  {DV_holo[0]:.3f}  (pull = {(DV_holo[0]-DV_bgs_data)/DV_bgs_err:+.2f} sigma)")
    print(f"  DESI w0waCDM: {DV_desi[0]:.3f}  (pull = {(DV_desi[0]-DV_bgs_data)/DV_bgs_err:+.2f} sigma)")


if __name__ == '__main__':
    create_combined_figure()