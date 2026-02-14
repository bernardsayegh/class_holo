import numpy as np
import matplotlib.pyplot as plt

# --- 1. Define the Models ---
H0_planck = 67.4
Om0_planck = 0.315
OL0_planck = 1.0 - Om0_planck

Om0_holo = 0.334
OL0_holo = 1.0 - Om0_holo
H0_holo = 67.4

def get_OL_z(z, Om0, OL0):
    E2 = Om0 * (1 + z)**3 + OL0
    return OL0 / E2

# --- 2. Real Observational Data ---
data_points = [
    (0.38, 81.2, 2.4),
    (0.51, 90.9, 2.3),
    (0.61, 99.0, 2.5),
    (0.698, 112.6, 2.2),
    (1.48, 152.0, 24.0),
    (2.33, 227.0, 8.0),
    (0.44, 82.6, 7.8),
    (0.60, 87.9, 6.1),
    (0.73, 97.3, 7.0),
    (0.0, 67.4, 0.5)
]

z_data = []
OL_data = []
yerr_data = []

for z_val, H_val, H_err in data_points:
    OL_val = OL0_planck * (H0_planck / H_val)**2
    OL_err = 2 * OL_val * (H_err / H_val)
    z_data.append(z_val)
    OL_data.append(OL_val)
    yerr_data.append(OL_err)

# --- 3. Plotting ---
z_range = np.linspace(0, 2.5, 500)
OL_holo_curve = get_OL_z(z_range, Om0_holo, OL0_holo)
OL_std_curve = get_OL_z(z_range, Om0_planck, OL0_planck)

fig, ax = plt.subplots(figsize=(12, 8))

ax.plot(z_range, OL_holo_curve, color='blue', linewidth=3.0, label=r'Holographic Model ($\Omega_m=0.334$)')
ax.plot(z_range, OL_std_curve, color='black', linestyle=':', linewidth=2.5, label=r'Standard $\Lambda$CDM (Planck)')

ax.axhline(y=1/3, color='black', linestyle='--', alpha=0.5, linewidth=2, label=r'Onset ($q=0$)')
ax.axhline(y=2/3, color='red', linestyle='--', alpha=0.5, linewidth=2, label='Holographic Saturation')

ax.errorbar(z_data, OL_data, yerr=yerr_data, fmt='o', color='orange', ecolor='orange', 
            markersize=9, capsize=5, elinewidth=2, markeredgewidth=2, label='Observational Data (BOSS/eBOSS/WiggleZ)')

ax.text(0.38, 0.52, r'BOSS $z=0.38$', color='orange', fontsize=14, ha='right', weight='bold')
ax.text(1.48, 0.15, r'eBOSS QSO', color='orange', fontsize=14, ha='left', weight='bold')
ax.text(2.33, 0.05, r'Ly-$\alpha$', color='orange', fontsize=14, ha='left', weight='bold')

ax.set_xlim(0, 2.5)
ax.set_ylim(0, 0.8)
ax.invert_xaxis()

ax.set_xlabel('Redshift $z$', fontsize=18)
ax.set_ylabel(r'Dark Energy Fraction $\Omega_{\Lambda}(z)$', fontsize=18)
ax.set_title('Holographic Bracketing with Expanded Observational Data', fontsize=20, pad=20)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True, linestyle=':', alpha=0.6, linewidth=1.5)
ax.legend(fontsize=14, loc='upper left', framealpha=0.9)

plt.tight_layout()
plt.savefig('holographic_data_bracket_expanded.png', dpi=300)
plt.savefig('holographic_data_bracket_expanded.pdf', dpi=300)
print("Saved: holographic_data_bracket_expanded.png/pdf")
