import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples
import os

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plt.rcParams.update({'font.size': 12})

path = "chains_schw/modelA_beta_free"
s = loadMCSamples(path, settings={"ignore_rows": 0.3})
param_names = [p.name for p in s.paramNames.names]
beta = s.samples[:, param_names.index("interaction_beta")]
sigma8 = s.samples[:, param_names.index("sigma8")]
Om = s.samples[:, param_names.index("Omega_m")]
S8 = sigma8 * np.sqrt(Om / 0.3)
S8_DES, S8_DES_ERR, S8_LCDM, BETA_THEORY = 0.776, 0.017, 0.834, 1.0/12.0
beta_mean = np.mean(beta)
idx = np.random.choice(len(beta), min(2000, len(beta)), replace=False)
ax1.scatter(beta[idx], S8[idx], c='blue', alpha=0.4, s=20, label='MCMC samples')
ax1.axhspan(S8_DES - S8_DES_ERR, S8_DES + S8_DES_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax1.axhline(S8_DES, color='green', ls='--', lw=2, label=f'DES Y3 $S_8$ = {S8_DES}')
ax1.axhline(S8_LCDM, color='red', ls=':', lw=2, label=f'LCDM $S_8$ = {S8_LCDM}')
ax1.axvline(BETA_THEORY, color='magenta', ls='--', lw=2, label=r'$\beta = 1/12$')
ax1.axvline(beta_mean, color='blue', ls=':', lw=1.5, label=f'MCMC mean beta = {beta_mean:.3f}')
ax1.set_xlabel(r'$\beta$', fontsize=14)
ax1.set_ylabel(r'$S_8$', fontsize=14)
ax1.set_xlim(0.03, 0.12)
ax1.set_ylim(0.74, 0.82)
ax1.set_title(r'(a) Entropic Brake: $\beta$--$S_8$ Anti-correlation', fontsize=13)
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, alpha=0.3)

configs = [("chains_schw/lcdm_shoes", "LCDM", "gray", "H0"), ("chains_schw/modelA_fixed", "Model A", "blue", "H0"), ("chains_schw/modelB_Amap2", "Model B", "red", "H0_local")]
samples = {}
for path, label, color, h0_param in configs:
    s = loadMCSamples(path, settings={"ignore_rows": 0.3})
    param_names = [p.name for p in s.paramNames.names]
    H0 = s.samples[:, param_names.index("H0_local")] if (h0_param == "H0_local" and "H0_local" in param_names) else s.samples[:, param_names.index("H0")]
    sigma8 = s.samples[:, param_names.index("sigma8")]
    Om = s.samples[:, param_names.index("Omega_m")]
    S8 = sigma8 * np.sqrt(Om / 0.3)
    samples[label] = {"H0": H0, "S8": S8, "color": color}

H0_SHOES, H0_SHOES_ERR = 73.04, 1.04
ax2.axvspan(H0_SHOES - H0_SHOES_ERR, H0_SHOES + H0_SHOES_ERR, color='khaki', alpha=0.4, zorder=0)
ax2.axvline(H0_SHOES, color='goldenrod', ls='--', lw=2, label=f'SH0ES $H_0$ = {H0_SHOES}')
ax2.axhspan(S8_DES - S8_DES_ERR, S8_DES + S8_DES_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax2.axhline(S8_DES, color='green', ls='--', lw=2, label=f'DES Y3 $S_8$ = {S8_DES}')
ax2.axhline(S8_LCDM, color='gray', ls=':', lw=1.5, alpha=0.7, label=f'LCDM $S_8$ = {S8_LCDM}')
for label, data in samples.items():
    idx = np.random.choice(len(data["H0"]), min(1500, len(data["H0"])), replace=False)
    ax2.scatter(data["H0"][idx], data["S8"][idx], c=data["color"], alpha=0.5, s=15, label=label)

ax2.set_xlabel(r'$H_0^{\rm eff}$ [km/s/Mpc]', fontsize=14)
ax2.set_ylabel(r'$S_8$', fontsize=14)
ax2.set_xlim(66, 75)
ax2.set_ylim(0.74, 0.85)
ax2.set_title(r'(b) $H_0$--$S_8$ Plane: Simultaneous Tension Resolution', fontsize=13)
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.annotate('LCDM\n(both tensions)', xy=(67.8, 0.82), fontsize=9, ha='center', color='gray', style='italic')
ax2.annotate('Model A\n(S8 resolved)', xy=(67.8, 0.775), fontsize=9, ha='center', color='blue', style='italic')
ax2.annotate('Model B\n(both resolved)', xy=(73.3, 0.775), fontsize=9, ha='center', color='red', style='italic')
plt.tight_layout()
plt.savefig("contours_combined.png", dpi=200, bbox_inches='tight')
plt.savefig("contours_combined.pdf", dpi=200, bbox_inches='tight')
print("Saved: contours_combined.png/pdf")
