#!/usr/bin/env python3
"""S8 vs clustering fraction with the v5 posterior overlay (fig:S8_panels).
Run from ~/class_holo_test (needs the CLASS python build on path)."""
import sys, os, glob; sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

BESTFIT = {
    'output': 'mPk', 'P_k_max_1/Mpc': 10.0,
    'omega_b': 0.02237, 'omega_cdm': 0.12, 'h': 0.68,
    'A_s': 2.0989032e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
    'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
}
lcdm = Class(); lcdm.set(BESTFIT); lcdm.compute()
S8_lcdm = lcdm.sigma8() * np.sqrt(lcdm.Omega_m() / 0.3)
lcdm.struct_cleanup(); lcdm.empty()

FCLUST_GRID = np.linspace(0, 1, 21)
S8_vals = []
for fc in FCLUST_GRID:
    params = dict(BESTFIT)
    params['interaction_beta'] = 1.0/12.0
    params['interaction_ieff_type'] = 4
    params['f_clust'] = fc
    c = Class(); c.set(params); c.compute()
    S8_vals.append(c.sigma8() * np.sqrt(c.Omega_m() / 0.3))
    c.struct_cleanup(); c.empty()
S8_vals = np.array(S8_vals)

# v5 posterior overlay: A(fc free) at beta=1/12
fc_mean = fc_std = None
CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")
cands = ["modelA_fcfree_v5_acc", "modelA_fcfree_v5"]
root = next((os.path.join(CHAINS_DIR, c) for c in cands
             if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt"))), None)
if root is None:
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, "modelA_fcfree*v5*.1.txt")))
    root = hits[0][:-6] if hits else None
if root:
    from getdist import loadMCSamples
    s = loadMCSamples(root, settings={"ignore_rows": 0.3})
    pn = [p.name for p in s.paramNames.names]
    key = "f_clust" if "f_clust" in pn else ("fclust" if "fclust" in pn else None)
    if key:
        w = s.weights; fcs = s.samples[:, pn.index(key)]
        fc_mean = np.average(fcs, weights=w)
        fc_std = np.sqrt(np.average((fcs - fc_mean)**2, weights=w))
        print(f"posterior overlay from {os.path.basename(root)}: "
              f"f_clust = {fc_mean:.3f} +/- {fc_std:.3f}  (grid: 0.20 +/- 0.16)")
else:
    print("note: no fcfree v5 chain found; plotting theory line only")

fig, ax = plt.subplots(figsize=(8, 6))
ax.axhspan(0.776 - 0.017, 0.776 + 0.017, alpha=0.15, color="green")
ax.axhline(0.776, color="green", ls="--", lw=2, label=r"DES Y3 $S_8 = 0.776 \pm 0.017$")
ax.axhline(S8_lcdm, color="gray", ls="--", lw=1.5, label=rf"CLASS $\Lambda$CDM $S_8 = {S8_lcdm:.3f}$")
ax.axhline(0.787, color="red", ls=":", lw=1.5, label=r"$\Lambda$CDM (joint fit) $S_8 = 0.787$")
if fc_mean is not None:
    ax.axvspan(fc_mean - fc_std, fc_mean + fc_std, color="blue", alpha=0.10, zorder=0)
    ax.axvline(fc_mean, color="blue", ls=":", lw=1.5,
               label=rf"MCMC $f_{{\rm clust}} = {fc_mean:.2f} \pm {fc_std:.2f}$")
ax.plot(FCLUST_GRID, S8_vals, lw=2.5, color="blue", label=r"Holographic $\beta = 1/12$")
ax.set_title(r"$S_8$ vs Clustering Fraction", fontsize=14)
ax.set_xlabel(r"$f_{\rm clust}$", fontsize=13)
ax.set_ylabel(r"$S_8 = \sigma_8\sqrt{\Omega_m/0.3}$", fontsize=13)
ax.set_xlim(-0.05, 1.05); ax.set_ylim(0.74, 0.86)
ax.legend(loc="upper left", fontsize=10)
ax.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig("S8_fclust.png", dpi=200, bbox_inches='tight')
plt.savefig("S8_fclust.pdf", bbox_inches='tight')
print(f"Saved: {os.path.abspath('S8_fclust.pdf')} (+ .png)")
print("\n=== caption batch values ===")
print(f"S8(fc=0) = {S8_vals[0]:.3f},  S8(fc=1) = {S8_vals[-1]:.3f},  CLASS LCDM = {S8_lcdm:.3f}")
