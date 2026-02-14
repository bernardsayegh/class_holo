import sys
sys.path.insert(0, 'python')
from classy import Class
from scipy.optimize import brentq
import numpy as np
import matplotlib.pyplot as plt

print("Creating combined Evolution + S(z) figure...")

THETA_S_TARGET = 1.040423
BETA_HOLO = 1.0/12.0

BASE_PARAMS = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
    "tau_reio": 0.054,
    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "Omega_k": 0.0,
    "lensing": "no",
}

def get_class_instance(h, beta):
    params = dict(BASE_PARAMS)
    params["h"] = float(h)
    if beta > 0:
        params["interaction_beta"] = float(beta)
        params["f_clust"] = 0.0
    c = Class()
    c.set(params)
    c.compute()
    return c

def solve_h(beta):
    def f(h):
        c = get_class_instance(h, beta)
        theta = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
        c.struct_cleanup()
        c.empty()
        return theta - THETA_S_TARGET
    return brentq(f, 0.55, 0.80, xtol=1e-6)

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

# Solve h
print("Solving h for LCDM...")
h_lcdm = solve_h(0.0)
print("Solving h for Holographic...")
h_holo = solve_h(BETA_HOLO)

# Get CLASS instances
c_lcdm = get_class_instance(h_lcdm, 0.0)
c_holo = get_class_instance(h_holo, BETA_HOLO)

Om0_l, Om0_h = c_lcdm.Omega_m(), c_holo.Omega_m()
Ode0_l, Ode0_h = 1.0 - Om0_l, 1.0 - Om0_h

bg_l = c_lcdm.get_background()
bg_h = c_holo.get_background()
c_lcdm.struct_cleanup(); c_lcdm.empty()
c_holo.struct_cleanup(); c_holo.empty()

z_l, Om_l, Ode_l = extract_omegas(bg_l)
z_h, Om_h, Ode_h = extract_omegas(bg_h)

# Mask to z <= 5
mask_l = z_l <= 5
mask_h = z_h <= 5
z_l, Om_l, Ode_l = z_l[mask_l], Om_l[mask_l], Ode_l[mask_l]
z_h, Om_h, Ode_h = z_h[mask_h], Om_h[mask_h], Ode_h[mask_h]

# S(z) = (9/2) * Omega_m * Omega_de
S_l = 4.5 * Om_l * Ode_l
S_h = 4.5 * Om_h * Ode_h

# Exact crossings from our calculation
z_enter = 0.609
z_exit = 0.011

print(f"LCDM: Om0={Om0_l:.3f}, Ode0={Ode0_l:.3f}")
print(f"Holo: Om0={Om0_h:.3f}, Ode0={Ode0_h:.3f}")

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# ============ Panel (a): Omega evolution ============
ax1.plot(z_l, Om_l, 'k--', lw=2.5, label=r'$\Lambda$CDM: $\Omega_m$')
ax1.plot(z_l, Ode_l, 'k--', lw=2.5, alpha=0.65, label=r'$\Lambda$CDM: $\Omega_{\rm de}$')
ax1.plot(z_h, Om_h, 'b-', lw=3, label=r'Holographic: $\Omega_m$')
ax1.plot(z_h, Ode_h, 'C1-', lw=3, label=r'Holographic: $\Omega_{\rm de}$')

# Horizontal guides
ax1.axhline(1/3, color='gray', ls='--', lw=1.8, alpha=0.8, label=r'$\Omega = 1/3$')
ax1.axhline(2/3, color='gray', ls='--', lw=1.8, alpha=0.8, label=r'$\Omega = 2/3$')

# Vertical lines
ax1.axvline(0.0, color='green', ls='-', lw=2.2, alpha=0.9, label=r'Present ($z=0$)')
ax1.axvline(0.61, color='purple', ls=':', lw=2.5, alpha=0.95, label=r'Holo: $\Omega_m=2/3$, $\Omega_{\rm de}=1/3$ at $z\approx 0.61$')

# Right-side annotations
ax1.text(0.12, Ode0_l + 0.08, rf'$\Lambda$CDM: $\Omega_{{\rm de}}(0)={Ode0_l:.3f}$', fontsize=11)
ax1.text(0.12, Ode0_h + 0.02, rf'Holo: $\Omega_{{\rm de}}(0)={Ode0_h:.3f}$', fontsize=11)
ax1.text(0.12, Om0_h - 0.05, rf'Holo: $\Omega_m(0)={Om0_h:.3f}$', fontsize=11)
ax1.text(0.12, Om0_l - 0.12, rf'$\Lambda$CDM: $\Omega_m(0)={Om0_l:.3f}$', fontsize=11)

ax1.set_xlabel('Redshift $z$', fontsize=14)
ax1.set_ylabel(r'Density fraction $\Omega$', fontsize=14)
ax1.set_title(r'(a) Density Evolution (matched to $100\theta_s$)', fontsize=13)
ax1.set_xlim(5, 0)
ax1.set_ylim(0, 1.06)
ax1.legend(loc='center', bbox_to_anchor=(0.52, 0.48), fontsize=9, framealpha=0.95)

# ============ Panel (b): S(z) window ============
# Mask for zoomed view
mask_l2 = z_l <= 2
mask_h2 = z_h <= 2
z_l2, S_l2 = z_l[mask_l2], S_l[mask_l2]
z_h2, S_h2 = z_h[mask_h2], S_h[mask_h2]

ax2.plot(z_l2, S_l2, 'k--', lw=2, label=r'$\Lambda$CDM')
ax2.plot(z_h2, S_h2, 'b-', lw=2.5, label=r'Holographic $\beta=1/12$')
ax2.axhline(1.0, color='red', ls='--', lw=2, label=r'$S = 1$ (Schwarzschild limit)')
ax2.fill_between(z_h2, S_h2, 1.0, where=(S_h2 >= 1), color='red', alpha=0.15, label='Super-Schwarzschild window')

# Vertical crossing lines
ax2.axvline(z_enter, color='purple', ls=':', lw=2, alpha=0.8)
ax2.axvline(z_exit, color='green', ls=':', lw=2, alpha=0.8)
ax2.text(z_enter - 0.05, 0.85, f'$z={z_enter:.2f}$', fontsize=11, color='purple', ha='right')
ax2.text(z_exit + 0.05, 0.85, f'$z={z_exit:.3f}$', fontsize=11, color='green', ha='left')

ax2.set_xlabel('Redshift $z$', fontsize=14)
ax2.set_ylabel(r'$S(z) = \frac{9}{2}\Omega_m\Omega_{\rm de}$', fontsize=14)
ax2.set_title('(b) Super-Schwarzschild Diagnostic', fontsize=13)
ax2.set_xlim(2, -0.05)
ax2.set_ylim(0.2, 1.15)
ax2.legend(loc='lower left', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('evolution_combined.png', dpi=200, bbox_inches='tight')
plt.savefig('evolution_combined.pdf', dpi=200, bbox_inches='tight')
print("Saved: evolution_combined.png/pdf")
