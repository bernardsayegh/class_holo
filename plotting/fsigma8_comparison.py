import sys; sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

rsd_data = [
    (0.067, 0.423, 0.055, 0.055, "6dFGS", "o"),
    (0.15,  0.490, 0.145, 0.145, "SDSS MGS", "s"),
    (0.38,  0.497, 0.045, 0.045, "BOSS DR12", "D"),
    (0.51,  0.458, 0.038, 0.038, "BOSS DR12", "D"),
    (0.61,  0.436, 0.034, 0.034, "BOSS DR12", "D"),
    (0.70,  0.473, 0.041, 0.041, "eBOSS LRG", "^"),
    (0.85,  0.315, 0.095, 0.095, "eBOSS ELG", "v"),
    (1.48,  0.462, 0.045, 0.045, "eBOSS QSO", "p"),
    (0.60,  0.550, 0.120, 0.120, "VIPERS", "*"),
    (0.86,  0.400, 0.110, 0.110, "VIPERS", "*"),
    (1.36,  0.482, 0.116, 0.116, "FastSound", "h"),
]

def compute_fsigma8(cosmo_params, z_arr, gamma=0.55, label=""):
    cosmo = Class()
    cosmo.set(cosmo_params)
    cosmo.compute()
    h = cosmo.h()
    R8 = 8.0 / h
    Om0 = cosmo.Omega_m()
    sig8_arr = np.array([cosmo.sigma(R8, z) for z in z_arr])

    # Extract actual Omega_m(z) from CLASS background
    bg = cosmo.get_background()
    z_bg = np.array(bg['z'], dtype=float)
    rho_crit = np.array(bg['(.)rho_crit'], dtype=float)
    rho_b = np.array(bg.get('(.)rho_b', np.zeros_like(z_bg)), dtype=float)
    rho_cdm = np.array(bg.get('(.)rho_cdm', np.zeros_like(z_bg)), dtype=float)
    Om_bg = (rho_b + rho_cdm) / rho_crit
    # Sort ascending for interp
    order = np.argsort(z_bg)
    Om_z = np.interp(z_arr, z_bg[order], Om_bg[order])

    f_arr = Om_z**gamma
    fsig8_arr = f_arr * sig8_arr

    print(f"\n{label}")
    print(f"  h = {h:.4f}, Om0 = {Om0:.4f}, sig8(0) = {sig8_arr[0]:.4f}")
    print(f"  {'z':>5}  {'sig8':>8}  {'Om(z)':>8}  {'f':>8}  {'fsig8':>8}")
    print("  " + "-" * 45)
    for zt in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
        s = np.interp(zt, z_arr, sig8_arr)
        o = np.interp(zt, z_arr, Om_z)
        fv = np.interp(zt, z_arr, f_arr)
        fs = np.interp(zt, z_arr, fsig8_arr)
        print(f"  {zt:5.2f}  {s:8.4f}  {o:8.4f}  {fv:8.4f}  {fs:8.4f}")

    cosmo.struct_cleanup()
    cosmo.empty()
    return fsig8_arr, sig8_arr, f_arr, Om_z

z_arr = np.linspace(0.02, 2.5, 300)

lcdm_params = {
    'omega_b': 0.02237, 'omega_cdm': 0.12, 'h': 0.68,
    'A_s': 2.0989032e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
    'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
    'output': 'mPk,tCl', 'P_k_max_1/Mpc': 10.0, 'z_max_pk': 3.5,
}

holo_params = dict(lcdm_params)
holo_params.update({
    'interaction_beta': 1.0/12.0,
    'interaction_ieff_type': 4,
    'f_clust': 0.0,
})

fsig8_l, sig8_l, f_l, Om_z_l = compute_fsigma8(lcdm_params, z_arr, gamma=0.55, label="LCDM")
fsig8_h, sig8_h, f_h, Om_z_h = compute_fsigma8(holo_params, z_arr, gamma=0.55, label="Holographic")

print("\n  Ratio (Holo / LCDM):")
print(f"  {'z':>5}  {'ratio':>10}  {'suppr':>12}")
print("  " + "-" * 30)
for zt in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
    fl = np.interp(zt, z_arr, fsig8_l)
    fh = np.interp(zt, z_arr, fsig8_h)
    print(f"  {zt:5.2f}  {fh/fl:10.4f}  {(1-fh/fl)*100:11.1f}%")

# Plot
c_lcdm = '#2166ac'
c_holo = '#b2182b'
c_scr  = '#fee08b'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7),
                                gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05},
                                sharex=True)

ax1.axvspan(0.02, 0.63, alpha=0.15, color=c_scr, label=r'SCR window ($\mathcal{S}>1$)')
ax2.axvspan(0.02, 0.63, alpha=0.15, color=c_scr)

ax1.plot(z_arr, fsig8_l, color=c_lcdm, lw=2.0, label=r'$\Lambda$CDM')
ax1.plot(z_arr, fsig8_h, color=c_holo, lw=2.0, ls='--', label=r'Holographic $\beta=1/12$')

surveys = {}
for z, fs8, eu, el, lab, mk in rsd_data:
    if lab not in surveys:
        surveys[lab] = {'z': [], 'fs8': [], 'eu': [], 'el': [], 'mk': mk}
    surveys[lab]['z'].append(z)
    surveys[lab]['fs8'].append(fs8)
    surveys[lab]['eu'].append(eu)
    surveys[lab]['el'].append(el)

colors_data = {
    '6dFGS': '#4daf4a', 'SDSS MGS': '#984ea3', 'BOSS DR12': '#ff7f00',
    'eBOSS LRG': '#a65628', 'eBOSS ELG': '#f781bf', 'eBOSS QSO': '#999999',
    'VIPERS': '#377eb8', 'FastSound': '#e41a1c',
}
for lab, d in surveys.items():
    c = colors_data.get(lab, 'gray')
    ax1.errorbar(d['z'], d['fs8'], yerr=[d['el'], d['eu']],
                 fmt=d['mk'], color=c, ms=7, capsize=3, lw=1.2, label=lab, zorder=5)

ax1.set_ylabel(r'$f\sigma_8(z)$', fontsize=14)
ax1.set_ylim(0.2, 0.65)
ax1.legend(fontsize=8, ncol=2, loc='upper right', framealpha=0.9)
ax1.tick_params(labelbottom=False)

ratio = fsig8_h / fsig8_l
ax2.plot(z_arr, ratio, color=c_holo, lw=2.0, ls='--')
ax2.axhline(1.0, color=c_lcdm, lw=1.0, ls=':')
ax2.fill_between([0.0, 2.2], 0.98, 1.02, alpha=0.1, color='gray',
                 label=r'Euclid/DESI $1\!-\!2\%$ forecast')

ax2.set_xlabel(r'Redshift $z$', fontsize=14)
ax2.set_ylabel(r'Holo / $\Lambda$CDM', fontsize=14)
ax2.set_ylim(0.94, 1.02)
ax2.set_xlim(0.0, 2.2)
ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.legend(fontsize=9, loc='lower right')
ax2.annotate(r'$\sim\!5\%$ suppression', xy=(0.08, 0.952),
            fontsize=9, color=c_holo, ha='left',
            arrowprops=dict(arrowstyle='->', color=c_holo, lw=1.2),
            xytext=(0.35, 0.943))

plt.tight_layout()
plt.savefig('fsigma8_comparison.png', dpi=200, bbox_inches='tight')
plt.savefig('fsigma8_comparison.pdf', dpi=300, bbox_inches='tight')
print("\nSaved: fsigma8_comparison.png/pdf")