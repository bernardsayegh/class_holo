#!/usr/bin/env python3
"""Regenerate posteriors.pdf (fig:posteriors) from the v5 publication chains.

Run from anywhere:  python3 posteriors_v5.py
Chains dir, root names, and the R25 anchor are v5-correct; the script
auto-resolves each chain root among known naming candidates and prints
posterior means +/- sigma for the caption batch.
"""
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples
from scipy.stats import gaussian_kde

CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")

# (candidate roots in priority order, label, color, which H0 to plot)
CONFIGS = [
    (["lcdm_v5_acc", "lcdm_v5_plancknu", "lcdm_v5"],
     "LCDM", "gray", "H0"),
    (["modelA_v5_acc", "modelA_v5"],
     r"Model A ($\beta=1/12$)", "blue", "H0"),
    (["modelB_A2_v5_acc", "modelB_A2_v5"],
     r"Model B ($A_{\rm map}=2$)", "red", "H0_local"),
]

def resolve_root(candidates):
    for c in candidates:
        if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt")):
            return os.path.join(CHAINS_DIR, c)
    # last resort: glob for anything matching the model prefix + v5
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, candidates[0].split("_v5")[0] + "*v5*.1.txt")))
    return hits[0][:-6] if hits else None

def col(s, *keys):
    pn = [p.name for p in s.paramNames.names]
    for k in keys:
        if k in pn:
            if k != keys[0]:
                print(f"  note: '{keys[0]}' not in chain; falling back to '{k}'")
            return s.samples[:, pn.index(k)]
    raise KeyError(f"none of {keys} in chain (have: {pn})")

samples = {}
for candidates, label, color, h0_key in CONFIGS:
    root = resolve_root(candidates)
    if root is None:
        print(f"WARNING: no v5 chain found for {candidates[0]} -- skipping {label}")
        continue
    s = loadMCSamples(root, settings={"ignore_rows": 0.3})
    H0     = col(s, h0_key, "H0")   # for Model B: warns loudly if H0_local is absent
    sigma8 = col(s, "sigma8")
    Om     = col(s, "Omega_m", "omegam")
    S8     = sigma8 * np.sqrt(Om / 0.3)
    samples[label] = dict(H0=H0, S8=S8, Omega_m=Om, sigma8=sigma8,
                          color=color, h0_key=h0_key)
    print(f"Loaded {label}: {os.path.basename(root)}  N={len(H0)}  (H0 column: {h0_key})")

# --- Observational targets (v5 conventions) ---
H0_SHOES, H0_SHOES_ERR = 73.49, 0.93        # Riess et al. 2025 (R25) -- was R22 73.04 +/- 1.04
S8_DES, S8_DES_ERR = 0.776, 0.017           # DES-Y3
OM_PLANCK, OM_PLANCK_ERR = 0.315, 0.007
OM_PANTHEON, OM_PANTHEON_ERR = 0.334, 0.018
SIGMA8_WL = S8_DES / np.sqrt(OM_PLANCK / 0.3)
SIGMA8_WL_ERR = S8_DES_ERR / np.sqrt(OM_PLANCK / 0.3)

def plot_kde(ax, data, color, label, xlim):
    x = np.linspace(*xlim, 500)
    y = gaussian_kde(data, bw_method=0.15)(x)
    ax.plot(x, y, color=color, lw=2.5, label=label)
    ax.fill_between(x, y, alpha=0.15, color=color)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

panels = [
    (axes[0, 0], "H0", (66, 76), r"$H_0^{\rm eff}$ [km/s/Mpc]",
     [(H0_SHOES, H0_SHOES_ERR, "lightgreen", "green", "-", "SH0ES (R25)")]),
    (axes[0, 1], "S8", (0.72, 0.86), r"$S_8$",
     [(S8_DES, S8_DES_ERR, "lightgreen", "green", "-", "DES-Y3")]),
    (axes[1, 0], "Omega_m", (0.28, 0.36), r"$\Omega_m$",
     [(OM_PLANCK, OM_PLANCK_ERR, "lightgreen", "green", "-", "Planck"),
      (OM_PANTHEON, OM_PANTHEON_ERR, "orange", "darkorange", "--", "Pantheon+ (shape)")]),
    (axes[1, 1], "sigma8", (0.72, 0.84), r"$\sigma_8$",
     [(SIGMA8_WL, SIGMA8_WL_ERR, "lightgreen", "green", "-", "DES-Y3")]),
]

for ax, key, xlim, xlabel, bands in panels:
    for c0, err, span_c, line_c, ls, blabel in bands:
        ax.axvspan(c0 - err, c0 + err, color=span_c, alpha=0.35, zorder=0)
        ax.axvline(c0, color=line_c, ls=ls, lw=2, label=blabel)
    for label, d in samples.items():
        plot_kde(ax, d[key], d["color"], label, xlim)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Probability density", fontsize=12)
    ax.set_xlim(xlim); ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="upper right")

plt.suptitle("Holographic Cosmology: Posterior Distributions (v5 grid)", fontsize=14, y=1.01)
plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"posteriors.{ext}", dpi=150, bbox_inches="tight")
print(f"\nSaved: {os.path.abspath('posteriors.pdf')} (+ .png)")

# --- Caption-batch summary: mean +/- sigma per model ---
print("\n=== posterior summary (for caption batch) ===")
for label, d in samples.items():
    line = "  ".join(f"{k}={np.mean(v):.3f}+/-{np.std(v):.3f}"
                     for k, v in d.items() if k in ("H0", "S8", "Omega_m", "sigma8"))
    print(f"{label} [{d['h0_key']}]: {line}")
plt.show()
