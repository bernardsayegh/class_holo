#!/usr/bin/env python3
"""S8/sigma8 barchart across models and weak-lensing surveys (fig:S8_barchart).
Run from anywhere: chains resolve from ~/class_holo_test/chains_test."""
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples

CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")
IGNORE_ROWS = 0.30

MODELS = [
    (["lcdm_v5_acc", "lcdm_v5_plancknu", "lcdm_v5"], r'$\Lambda$CDM'),
    (["modelA_v5_acc", "modelA_v5"], "Model A\n" + r'($\beta=1/12$)'),
    (["modelB_A2_v5_acc", "modelB_A2_v5"], "Model B\n" + r'($A_{\rm map}=2$)'),
]

OBS = {
    "DES Y3":         {"S8": 0.776, "S8_err": 0.017, "sigma8": 0.760, "sigma8_err": 0.020},
    "DES Y6\n(NLA)":  {"S8": 0.798, "S8_err": 0.015, "sigma8": None, "sigma8_err": None},
    "DES Y6\n(TATT)": {"S8": 0.783, "S8_err": 0.017, "sigma8": None, "sigma8_err": None},
    "KiDS-1000":      {"S8": 0.759, "S8_err": 0.024, "sigma8": None, "sigma8_err": None},
    "KiDS-Leg.":      {"S8": 0.815, "S8_err": 0.019, "sigma8": None, "sigma8_err": None},
    "HSC-Y3":         {"S8": 0.776, "S8_err": 0.033, "sigma8": None, "sigma8_err": None},
}

CLASS_THETAS = {
    r'$\Lambda$CDM': {"sigma8": 0.80761, "S8": 0.82841},
    "Model A\n" + r'($\beta=1/12$)': {"sigma8": 0.74255, "S8": 0.77793},
}

def find_chain_base(candidates):
    for c in candidates:
        if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt")):
            return os.path.join(CHAINS_DIR, c)
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, candidates[0].split("_v5")[0] + "*v5*.1.txt")))
    if hits:
        return hits[0][:-6]
    raise FileNotFoundError(f"no v5 chain for {candidates[0]} in {CHAINS_DIR}")

def ensure_derived_S8(s):
    names = [p.name for p in s.paramNames.names]
    if "S8" in names:
        return
    sig8 = s.samples[:, names.index("sigma8")]
    Om = s.samples[:, names.index("Omega_m")]
    s.addDerived(sig8 * np.sqrt(Om / 0.3), name="S8", label="S_8")

def get_stats(s, param):
    names = [p.name for p in s.paramNames.names]
    i = names.index(param)
    return float(s.getMeans()[i]), float(np.sqrt(s.getCov()[i, i]))

results = {}
for candidates, label in MODELS:
    base = find_chain_base(candidates)
    s = loadMCSamples(base, settings={"ignore_rows": IGNORE_ROWS})
    ensure_derived_S8(s)
    sig8, sig8_err = get_stats(s, "sigma8")
    S8, S8_err = get_stats(s, "S8")
    results[label] = dict(sigma8=sig8, sigma8_err=sig8_err, S8=S8, S8_err=S8_err)
    print(f"Loaded {label.splitlines()[0]}: {os.path.basename(base)}  "
          f"sigma8={sig8:.3f}+/-{sig8_err:.3f}  S8={S8:.3f}+/-{S8_err:.3f}")

model_labels = [m[1] for m in MODELS]
labels_a = model_labels + list(OBS.keys())
n = len(labels_a)

sig8_vals = [results[l]["sigma8"] for l in model_labels]
sig8_errs = [results[l]["sigma8_err"] for l in model_labels]
S8_vals   = [results[l]["S8"] for l in model_labels]
S8_errs   = [results[l]["S8_err"] for l in model_labels]
for k in OBS:
    sig8_vals.append(OBS[k]["sigma8"] if OBS[k]["sigma8"] else np.nan)
    sig8_errs.append(OBS[k]["sigma8_err"] if OBS[k]["sigma8_err"] else 0.0)
    S8_vals.append(OBS[k]["S8"])
    S8_errs.append(OBS[k]["S8_err"])

fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(n)
w = 0.30

ax.axhspan(0.776 - 0.017, 0.776 + 0.017, alpha=0.15, color="green",
           label=r"DES Y3 $S_8$ band")
ax.axhline(0.834, color='grey', ls='--', lw=1.0, alpha=0.6,
           label=r'Planck $\Lambda$CDM $S_8 = 0.834$')

bars_sig8 = ax.bar(x - w/2, sig8_vals, w, yerr=sig8_errs, capsize=4,
                   color='steelblue', edgecolor='steelblue', label=r"$\sigma_8$")
bars_S8 = ax.bar(x + w/2, S8_vals, w, yerr=S8_errs, capsize=4,
                 color='coral', edgecolor='coral', label=r"$S_8$")

for i, v in enumerate(sig8_vals):
    if np.isnan(v):
        bars_sig8[i].set_alpha(0.0)

for b, v, e in zip(bars_sig8, sig8_vals, sig8_errs):
    if not np.isnan(v):
        ax.text(b.get_x() + b.get_width()/2, v + e + 0.004,
                f"{v:.3f}", ha="center", va="bottom", fontsize=8.5)
for b, v, e in zip(bars_S8, S8_vals, S8_errs):
    ax.text(b.get_x() + b.get_width()/2, v + e + 0.004,
            f"{v:.3f}", ha="center", va="bottom", fontsize=8.5)

for i, lab in enumerate(model_labels):
    cv = CLASS_THETAS.get(lab)
    if cv:
        ax.plot(i - w/2, cv["sigma8"], 'ko', ms=7, zorder=6)
        ax.plot(i + w/2, cv["S8"], 'ko', ms=7, zorder=6)
ax.plot([], [], 'ko', ms=7, label=r"CLASS ($\theta_s$-matched)")

ax.axvline(x=len(MODELS) - 0.5, color='grey', ls=':', lw=0.8, alpha=0.5)
ax.set_title("Clustering Amplitude Comparison", fontsize=14)
ax.set_ylabel("Clustering Amplitude", fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(labels_a, fontsize=10)
ax.set_ylim(0.70, 0.88)
ax.legend(loc="upper right", fontsize=9, framealpha=0.95)

plt.tight_layout()
plt.savefig("S8_barchart.png", dpi=200, bbox_inches='tight')
plt.savefig("S8_barchart.pdf", bbox_inches='tight')
print(f"\nSaved: {os.path.abspath('S8_barchart.pdf')} (+ .png)")

print("\n=== caption batch values ===")
for lab in model_labels:
    r = results[lab]
    print(f"{lab.splitlines()[0]}: sigma8 = {r['sigma8']:.3f} +/- {r['sigma8_err']:.3f},  "
          f"S8 = {r['S8']:.3f} +/- {r['S8_err']:.3f}")
