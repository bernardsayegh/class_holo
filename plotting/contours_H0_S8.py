#!/usr/bin/env python3
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples

np.random.seed(42)
plt.rcParams.update({'font.size': 12})
S8_DES, S8_DES_ERR = 0.776, 0.017
S8_PLANCK_LCDM = 0.834
CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")

def find_root(candidates):
    for c in candidates:
        if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt")):
            return os.path.join(CHAINS_DIR, c)
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, candidates[0].split("_v5")[0] + "*v5*.1.txt")))
    if hits:
        return hits[0][:-6]
    raise FileNotFoundError(f"no v5 chain for {candidates[0]}")

configs = [
    (["lcdm_v5_acc", "lcdm_v5_plancknu", "lcdm_v5"], "LCDM", "gray", "H0"),
    (["modelA_v5_acc", "modelA_v5"], "Model A", "blue", "H0"),
    (["modelB_A2_v5_acc", "modelB_A2_v5"], "Model B", "red", "H0_local"),
]
samples = {}
for candidates, label, color, h0_param in configs:
    root = find_root(candidates)
    s = loadMCSamples(root, settings={"ignore_rows": 0.3})
    param_names = [p.name for p in s.paramNames.names]
    H0 = s.samples[:, param_names.index("H0_local")] if (h0_param == "H0_local" and "H0_local" in param_names) else s.samples[:, param_names.index("H0")]
    sigma8 = s.samples[:, param_names.index("sigma8")]
    Om = s.samples[:, param_names.index("Omega_m")]
    S8 = sigma8 * np.sqrt(Om / 0.3)
    w = s.weights
    samples[label] = {"H0": H0, "S8": S8, "color": color, "w": w}
    print(f"Loaded {label}: {os.path.basename(root)}  H0={np.average(H0,weights=w):.2f}  S8={np.average(S8,weights=w):.3f}")

H0_SHOES, H0_SHOES_ERR = 73.49, 0.93   # SH0ES R25

fig, ax = plt.subplots(figsize=(8, 6))
ax.axvspan(H0_SHOES - H0_SHOES_ERR, H0_SHOES + H0_SHOES_ERR, color='khaki', alpha=0.4, zorder=0)
ax.axvline(H0_SHOES, color='goldenrod', ls='--', lw=2, label=f'SH0ES (R25) $H_0$ = {H0_SHOES}')
ax.axhspan(S8_DES - S8_DES_ERR, S8_DES + S8_DES_ERR, color='lightgreen', alpha=0.4, zorder=0)
ax.axhline(S8_DES, color='green', ls='--', lw=2, label=f'DES Y3 $S_8$ = {S8_DES}')
ax.axhline(S8_PLANCK_LCDM, color='gray', ls=':', lw=1.5, alpha=0.7, label=fr'Planck $\Lambda$CDM $S_8$ = {S8_PLANCK_LCDM}')
for label, data in samples.items():
    prob = data["w"] / data["w"].sum()
    idx = np.random.choice(len(data["H0"]), min(1500, len(data["H0"])), replace=True, p=prob)
    ax.scatter(data["H0"][idx], data["S8"][idx], c=data["color"], alpha=0.5, s=15, label=label)
ax.set_xlabel(r'$H_0^{\rm eff}$ [km/s/Mpc]', fontsize=14)
ax.set_ylabel(r'$S_8$', fontsize=14)
ax.set_xlim(66, 76); ax.set_ylim(0.72, 0.86)
ax.set_title(r'$H_0$--$S_8$ Plane: Simultaneous Tension Resolution', fontsize=13)
ax.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)
ax.annotate('LCDM\n(both tensions)', xy=(69.1, 0.800), fontsize=9, ha='center', color='gray', style='italic')
ax.annotate('Model A\n($S_8$ resolved)', xy=(67.4, 0.752), fontsize=9, ha='center', color='blue', style='italic')
ax.annotate('Model B\n(both resolved)', xy=(73.6, 0.757), fontsize=9, ha='center', color='red', style='italic')
plt.tight_layout()
plt.savefig("contours_H0_S8.png", dpi=200, bbox_inches='tight')
plt.savefig("contours_H0_S8.pdf", dpi=200, bbox_inches='tight')
print("Saved: contours_H0_S8.png/pdf")
