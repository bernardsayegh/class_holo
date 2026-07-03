#!/usr/bin/env python3
"""Entropic-brake beta--S8 scatter, Model A beta-free (fig:S8_contours).
Weighted samples; run from anywhere."""
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples

np.random.seed(42)
plt.rcParams.update({'font.size': 12})
CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")

def find_root(candidates):
    for c in candidates:
        if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt")):
            return os.path.join(CHAINS_DIR, c)
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, candidates[0].split("_v5")[0] + "*v5*.1.txt")))
    if hits:
        return hits[0][:-6]
    raise FileNotFoundError(f"no v5 chain for {candidates[0]} in {CHAINS_DIR}")

def col(s, *keys):
    pn = [p.name for p in s.paramNames.names]
    for k in keys:
        if k in pn:
            return s.samples[:, pn.index(k)]
    raise KeyError(f"none of {keys} in chain (have: {pn})")

root = find_root(["modelA_bfree_v5_acc", "modelA_bfree_v5"])
s = loadMCSamples(root, settings={"ignore_rows": 0.3})
weights = s.weights
beta   = col(s, "interaction_beta", "beta")
sigma8 = col(s, "sigma8")
Om     = col(s, "Omega_m", "omegam")
S8 = sigma8 * np.sqrt(Om / 0.3)

beta_mean = np.average(beta, weights=weights)
beta_std  = np.sqrt(np.average((beta - beta_mean)**2, weights=weights))
S8_mean = np.average(S8, weights=weights)
S8_std  = np.sqrt(np.average((S8 - S8_mean)**2, weights=weights))
print(f"Loaded: {os.path.basename(root)}")
print(f"beta = {beta_mean:.4f} +/- {beta_std:.4f}   (grid: 0.102 +/- 0.043)")
print(f"S8   = {S8_mean:.4f} +/- {S8_std:.4f}")

S8_LCDM_JOINT = 0.787   # v5 anchored LCDM joint fit

fig, ax = plt.subplots(figsize=(9, 6))
prob = weights / weights.sum()
idx = np.random.choice(len(beta), min(2000, len(beta)), replace=True, p=prob)
ax.scatter(beta[idx], S8[idx], c='blue', alpha=0.4, s=20, label='MCMC samples')
ax.axhspan(0.776 - 0.017, 0.776 + 0.017, color='lightgreen', alpha=0.4, zorder=0)
ax.axhline(0.776, color='green', ls='--', lw=2, label='DES Y3 $S_8$ = 0.776')
ax.axhline(S8_LCDM_JOINT, color='red', ls=':', lw=2,
           label=fr'$\Lambda$CDM (joint fit) $S_8$ = {S8_LCDM_JOINT}')
ax.axvline(1.0/12.0, color='magenta', ls='--', lw=2, label=r'$\beta = 1/12$')
ax.axvline(beta_mean, color='blue', ls=':', lw=1.5,
           label='MCMC mean $\\beta$ = %.3f' % beta_mean)
ax.set_xlabel(r'$\beta$', fontsize=14)
ax.set_ylabel(r'$S_8$', fontsize=14)
xlo = max(0.0, np.percentile(beta, 0.5) - 0.01)
xhi = np.percentile(beta, 99.5) + 0.01
ylo = np.percentile(S8, 0.5) - 0.012
yhi = max(np.percentile(S8, 99.5), S8_LCDM_JOINT) + 0.012
ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
ax.set_title(r'Entropic Brake: $\beta$--$S_8$ Anti-correlation', fontsize=13)
ax.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("contours_beta_S8.png", dpi=200, bbox_inches='tight')
plt.savefig("contours_beta_S8.pdf", bbox_inches='tight')
print(f"\nSaved: {os.path.abspath('contours_beta_S8.pdf')} (+ .png)")
print("\n=== caption batch values ===")
print(f"beta = {beta_mean:.3f} +/- {beta_std:.3f};  LCDM joint-fit line = {S8_LCDM_JOINT}")
