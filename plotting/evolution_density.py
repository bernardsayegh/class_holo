import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

BESTFIT = {
    "output": "tCl",
    "omega_b": 0.02237, "omega_cdm": 0.12, "h": 0.68,
    "A_s": 2.0989032e-9, "n_s": 0.9649, "tau_reio": 0.0544,
    "N_ur": 2.0328, "N_ncdm": 1, "m_ncdm": 0.06,
    "Omega_k": 0.0, "YHe": 0.2454, "lensing": "no",
}

lcdm = Class(); lcdm.set(BESTFIT); lcdm.compute()
params_h = dict(BESTFIT)
params_h["interaction_beta"] = 1.0/12.0
params_h["interaction_ieff_type"] = 4
params_h["f_clust"] = 0.0
holo = Class(); holo.set(params_h); holo.compute()

Om0_l, Om0_h = lcdm.Omega_m(), holo.Omega_m()
Ode0_l, Ode0_h = 1.0 - Om0_l, 1.0 - Om0_h
print(f"LCDM:  Omega_m(0) = {Om0_l:.4f}, Omega_de(0) = {Ode0_l:.4f}")
print(f"Holo:  Omega_m(0) = {Om0_h:.4f}, Omega_de(0) = {Ode0_h:.4f}")

def extract_omegas(bg):
    z = np.array(bg["z"], dtype=float)
    rho_crit = np.array(bg["(.)rho_crit"], dtype=float)
    rho_b = np.array(bg.get("(.)rho_b", np.zeros_like(z)), dtype=float)
    rho_cdm = np.array(bg.get("(.)rho_cdm", np.zeros_like(z)), dtype=float)
    Om = (rho_b + rho_cdm) / rho_crit
    if "(.)rho_fld" in bg:
        Ode = np.array(bg["(.)rho_fld"], dtype=float) / rho_crit
    elif "(.)rho_lambda" in bg:
        Ode = np.array(bg["(.)rho_lambda"], dtype=float) / rho_crit
    else:
        Ode = 1.0 - Om
    order = np.argsort(z)
    return z[order], Om[order], Ode[order]

bg_l = lcdm.get_background(); bg_h = holo.get_background()
lcdm.struct_cleanup(); lcdm.empty()
holo.struct_cleanup(); holo.empty()

z_l, Om_l, Ode_l = extract_omegas(bg_l)
z_h, Om_h, Ode_h = extract_omegas(bg_h)

mask_l = z_l <= 5; mask_h = z_h <= 5
z_l, Om_l, Ode_l = z_l[mask_l], Om_l[mask_l], Ode_l[mask_l]
z_h, Om_h, Ode_h = z_h[mask_h], Om_h[mask_h], Ode_h[mask_h]

# Find SCR entry: S = 9/2 Om Ode = 1 (consistent with S(z) script)
S_h = 4.5 * Om_h * Ode_h
for i in range(len(z_h) - 1):
    if (S_h[i] - 1.0) * (S_h[i+1] - 1.0) < 0 and z_h[i] > 0.1:
        z_cross = z_h[i] + (1.0 - S_h[i]) * (z_h[i+1] - z_h[i]) / (S_h[i+1] - S_h[i])
        print(f"SCR entry (S=1) at z = {z_cross:.4f}")
        break

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(z_l, Om_l, 'k--', lw=2.5, label=r'$\Lambda$CDM: $\Omega_m$')
ax.plot(z_l, Ode_l, 'k--', lw=2.5, alpha=0.65, label=r'$\Lambda$CDM: $\Omega_{\rm de}$')
ax.plot(z_h, Om_h, 'b-', lw=3, label=r'Holographic: $\Omega_m$')
ax.plot(z_h, Ode_h, 'C1-', lw=3, label=r'Holographic: $\Omega_{\rm de}$')
ax.axhline(1/3, color='gray', ls='--', lw=1.8, alpha=0.8, label=r'$\Omega = 1/3$')
ax.axhline(2/3, color='gray', ls='--', lw=1.8, alpha=0.8, label=r'$\Omega = 2/3$')
ax.axvline(0.0, color='green', ls='-', lw=2.2, alpha=0.9, label=r'Present ($z=0$)')
ax.axvline(z_cross, color='purple', ls=':', lw=2.5, alpha=0.95,
           label=rf'SCR entry: $\mathcal{{S}}=1$ at $z\approx {z_cross:.2f}$')
ax.text(0.12, Ode0_l + 0.08, rf'$\Lambda$CDM: $\Omega_{{\rm de}}(0)={Ode0_l:.3f}$', fontsize=11)
ax.text(0.12, Ode0_h + 0.02, rf'Holo: $\Omega_{{\rm de}}(0)={Ode0_h:.3f}$', fontsize=11)
ax.text(0.12, Om0_h - 0.05, rf'Holo: $\Omega_m(0)={Om0_h:.3f}$', fontsize=11)
ax.text(0.12, Om0_l - 0.12, rf'$\Lambda$CDM: $\Omega_m(0)={Om0_l:.3f}$', fontsize=11)
ax.set_xlabel('Redshift $z$', fontsize=14)
ax.set_ylabel(r'Density fraction $\Omega$', fontsize=14)
ax.set_title(r'Density Evolution (Planck 2018 cosmology)', fontsize=13)
ax.set_xlim(5, 0); ax.set_ylim(0, 1.06)
ax.legend(loc='center', bbox_to_anchor=(0.52, 0.48), fontsize=9, framealpha=0.95)
plt.tight_layout()
plt.savefig('evolution_density.png', dpi=200, bbox_inches='tight')
plt.savefig('evolution_density.pdf', dpi=200, bbox_inches='tight')
print("\nSaved: evolution_density.png/pdf")