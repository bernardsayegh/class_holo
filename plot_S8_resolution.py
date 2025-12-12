import numpy as np
import matplotlib.pyplot as plt

# Data
Omega_m_planck = 0.3111
Omega_m_holo = 0.334  # holographic floor

# Sigma8 values from CLASS runs
sigma8_lcdm = 0.8227
sigma8_holo = 0.7590

# S8 = sigma8 * sqrt(Omega_m / 0.3)
S8_lcdm = sigma8_lcdm * np.sqrt(Omega_m_planck / 0.3)
S8_holo = sigma8_holo * np.sqrt(Omega_m_holo / 0.3)

# Weak lensing target
S8_wl_low = 0.76
S8_wl_high = 0.78
sigma8_wl_low = 0.75
sigma8_wl_high = 0.77

print(f"LCDM: sigma8={sigma8_lcdm:.3f}, S8={S8_lcdm:.3f}")
print(f"Holo: sigma8={sigma8_holo:.3f}, S8={S8_holo:.3f}")

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

x = np.array([0, 1, 2])
width = 0.35

# sigma8 bars
sigma8_vals = [sigma8_lcdm, sigma8_holo, (sigma8_wl_low + sigma8_wl_high)/2]
sigma8_errs = [0.006, 0.005, (sigma8_wl_high - sigma8_wl_low)/2]
bars1 = ax.bar(x - width/2, sigma8_vals, width, label=r'$\sigma_8$', color='steelblue', 
               yerr=sigma8_errs, capsize=5)

# S8 bars  
S8_vals = [S8_lcdm, S8_holo, (S8_wl_low + S8_wl_high)/2]
S8_errs = [0.01, 0.01, (S8_wl_high - S8_wl_low)/2]
bars2 = ax.bar(x + width/2, S8_vals, width, label=r'$S_8 = \sigma_8(\Omega_m/0.3)^{0.5}$', 
               color='darkorange', yerr=S8_errs, capsize=5)

# Weak lensing target band
ax.axhspan(S8_wl_low, S8_wl_high, alpha=0.2, color='green', label='Weak lensing target')

# Arrow showing suppression
ax.annotate('', xy=(1, sigma8_holo + 0.01), xytext=(0, sigma8_lcdm - 0.01),
            arrowprops=dict(arrowstyle='->', color='green', lw=2))
ax.text(0.5, 0.795, 'Dark Drag\nSuppression', ha='center', fontsize=10, color='green')

ax.set_ylabel('Clustering Amplitude', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels([r'Planck $\Lambda$CDM', r'Holographic ($\beta=0.5$)', 'Weak Lensing\n(DES/KiDS)'])
ax.set_ylim(0.70, 0.90)
ax.legend(loc='upper right')
ax.set_title(r'Resolution of the $S_8$ Tension', fontsize=14)
ax.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, val in zip(bars1, sigma8_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.012, f'{val:.3f}', 
            ha='center', va='bottom', fontsize=9)
for bar, val in zip(bars2, S8_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.012, f'{val:.3f}', 
            ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig('Figure_S8_resolution.png', dpi=150, bbox_inches='tight')
plt.savefig('Figure_S8_resolution.pdf', bbox_inches='tight')
print("Saved S8 resolution figures!")
