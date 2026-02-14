import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples
from scipy.stats import gaussian_kde
import os

# Load chains
configs = [
    ("chains_schw/lcdm_shoes", "LCDM", "gray", "H0"),
    ("chains_schw/modelA_fixed", "Model A (beta=1/12)", "blue", "H0"),
    ("chains_schw/modelB_Amap2", "Model B (Amap=2)", "red", "H0_local"),
]

samples = {}
for path, label, color, h0_param in configs:
    if os.path.exists(f"{path}.1.txt"):
        s = loadMCSamples(path, settings={"ignore_rows": 0.3})
        param_names = [p.name for p in s.paramNames.names]

        if h0_param == "H0_local" and "H0_local" in param_names:
            H0 = s.samples[:, param_names.index("H0_local")]
        else:
            H0 = s.samples[:, param_names.index("H0")]

        sigma8 = s.samples[:, param_names.index("sigma8")]
        Om = s.samples[:, param_names.index("Omega_m")]
        S8 = sigma8 * np.sqrt(Om / 0.3)

        samples[label] = {"H0": H0, "S8": S8, "Omega_m": Om, "sigma8": sigma8,
                          "color": color, "h0_param": h0_param}
        print(f"Loaded {label}: N={len(H0)}, using {h0_param}")
    else:
        print(f"WARNING: {path}.1.txt not found, skipping")

# Observational targets
H0_SHOES = 73.04
H0_SHOES_ERR = 1.04
S8_DES = 0.776
S8_DES_ERR = 0.017
OM_PLANCK = 0.315
OM_PLANCK_ERR = 0.007
OM_PANTHEON = 0.334
OM_PANTHEON_ERR = 0.018
SIGMA8_WL = S8_DES / np.sqrt(OM_PLANCK / 0.3)
SIGMA8_WL_ERR = S8_DES_ERR / np.sqrt(OM_PLANCK / 0.3)

def plot_kde(ax, data, color, label, xlim):
    x = np.linspace(xlim[0], xlim[1], 500)
    kde = gaussian_kde(data, bw_method=0.15)
    y = kde(x)
    ax.plot(x, y, color=color, lw=2.5, label=label)
    ax.fill_between(x, y, alpha=0.15, color=color)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# --- Panel 1: H0 ---
ax = axes[0, 0]
xlim = (66, 76)
ax.axvspan(H0_SHOES - H0_SHOES_ERR, H0_SHOES + H0_SHOES_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax.axvline(H0_SHOES, color='green', ls='-', lw=2, label='SH0ES')
for label, data in samples.items():
    plot_kde(ax, data["H0"], data["color"], label, xlim)
ax.set_xlabel(r"$H_0^{\rm eff}$ [km/s/Mpc]", fontsize=12)
ax.set_ylabel("Probability density", fontsize=12)
ax.set_xlim(xlim)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9, loc='upper right')

# --- Panel 2: S8 ---
ax = axes[0, 1]
xlim = (0.72, 0.86)
ax.axvspan(S8_DES - S8_DES_ERR, S8_DES + S8_DES_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax.axvline(S8_DES, color='green', ls='-', lw=2, label='DES-Y3')
for label, data in samples.items():
    plot_kde(ax, data["S8"], data["color"], label, xlim)
ax.set_xlabel(r"$S_8$", fontsize=12)
ax.set_ylabel("Probability density", fontsize=12)
ax.set_xlim(xlim)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9, loc='upper right')

# --- Panel 3: Omega_m ---
ax = axes[1, 0]
xlim = (0.28, 0.36)
ax.axvspan(OM_PLANCK - OM_PLANCK_ERR, OM_PLANCK + OM_PLANCK_ERR, color='lightgreen', alpha=0.3, zorder=0)
ax.axvline(OM_PLANCK, color='green', ls='-', lw=2, label='Planck')
ax.axvspan(OM_PANTHEON - OM_PANTHEON_ERR, OM_PANTHEON + OM_PANTHEON_ERR, color='orange', alpha=0.2, zorder=0)
ax.axvline(OM_PANTHEON, color='darkorange', ls='--', lw=2, label='Pantheon+ (shape)')
for label, data in samples.items():
    plot_kde(ax, data["Omega_m"], data["color"], label, xlim)
ax.set_xlabel(r"$\Omega_m$", fontsize=12)
ax.set_ylabel("Probability density", fontsize=12)
ax.set_xlim(xlim)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9, loc='upper right')

# --- Panel 4: sigma8 ---
ax = axes[1, 1]
xlim = (0.72, 0.84)
ax.axvspan(SIGMA8_WL - SIGMA8_WL_ERR, SIGMA8_WL + SIGMA8_WL_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax.axvline(SIGMA8_WL, color='green', ls='-', lw=2, label='DES-Y3')
for label, data in samples.items():
    plot_kde(ax, data["sigma8"], data["color"], label, xlim)
ax.set_xlabel(r"$\sigma_8$", fontsize=12)
ax.set_ylabel("Probability density", fontsize=12)
ax.set_xlim(xlim)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9, loc='upper right')

plt.suptitle("Holographic Cosmology: Posterior Distributions", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig("posteriors.png", dpi=150, bbox_inches='tight')
plt.savefig("posteriors.pdf", dpi=150, bbox_inches='tight')
print("\nSaved: posteriors.png and posteriors.pdf")
plt.show()
