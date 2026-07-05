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

def find_crossings(z, S, threshold=1.0):
    crossings = []
    for i in range(len(z) - 1):
        if (S[i] - threshold) * (S[i+1] - threshold) < 0:
            zc = z[i] + (threshold - S[i]) * (z[i+1] - z[i]) / (S[i+1] - S[i])
            crossings.append(zc)
    return sorted(crossings, reverse=True)

bg_l = lcdm.get_background(); bg_h = holo.get_background()
lcdm.struct_cleanup(); lcdm.empty()
holo.struct_cleanup(); holo.empty()

z_l, Om_l, Ode_l = extract_omegas(bg_l)
z_h, Om_h, Ode_h = extract_omegas(bg_h)

mask_l = z_l <= 2; mask_h = z_h <= 2
z_l, Om_l, Ode_l = z_l[mask_l], Om_l[mask_l], Ode_l[mask_l]
z_h, Om_h, Ode_h = z_h[mask_h], Om_h[mask_h], Ode_h[mask_h]

S_l = 4.5 * Om_l * Ode_l
S_h = 4.5 * Om_h * Ode_h

crossings_l = find_crossings(z_l, S_l)
crossings_h = find_crossings(z_h, S_h)

print(f"LCDM crossings at S=1:    {['z={:.4f}'.format(zc) for zc in crossings_l]}")
print(f"Holographic crossings at S=1: {['z={:.4f}'.format(zc) for zc in crossings_h]}")

z_enter = crossings_h[0] if len(crossings_h) >= 1 else 0.609
z_exit  = crossings_h[1] if len(crossings_h) >= 2 else 0.011
print(f"\nSCR window: z_enter = {z_enter:.4f}, z_exit = {z_exit:.4f}")
print(f"S_max(LCDM) = {S_l.max():.4f} at z = {z_l[np.argmax(S_l)]:.3f}")
print(f"S_max(holo)  = {S_h.max():.4f} at z = {z_h[np.argmax(S_h)]:.3f}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(z_l, S_l, 'k--', lw=2, label=r'$\Lambda$CDM')
ax.plot(z_h, S_h, 'b-', lw=2.5, label=r'Holographic $\beta=1/12$')
ax.axhline(1.0, color='red', ls='--', lw=2, label=r'$\mathcal{S} = 1$ (Schwarzschild limit)')
ax.fill_between(z_h, S_h, 1.0, where=(S_h >= 1), color='red', alpha=0.15,
                label='Super-Schwarzschild window')
ax.axvline(z_enter, color='purple', ls=':', lw=2, alpha=0.8)
ax.axvline(z_exit,  color='green',  ls=':', lw=2, alpha=0.8)
ax.text(z_enter - 0.05, 0.85, f'$z={z_enter:.2f}$', fontsize=11, color='purple', ha='right')
ax.text(z_exit + 0.05,  0.85, f'$z={z_exit:.3f}$',  fontsize=11, color='green',  ha='left')
ax.set_xlabel('Redshift $z$', fontsize=14)
ax.set_ylabel(r'$\mathcal{S}(z) = \frac{9}{2}\,\Omega_m\,\Omega_{\rm de}$', fontsize=14)
ax.set_title('Super-Schwarzschild Diagnostic (Planck 2018 cosmology)', fontsize=13)
ax.set_xlim(2, -0.05); ax.set_ylim(0.2, 1.15)
ax.legend(loc='lower left', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('evolution_Sz.png', dpi=200, bbox_inches='tight')
plt.savefig('evolution_Sz.pdf', dpi=200, bbox_inches='tight')
print("\nSaved: evolution_Sz.png/pdf")
