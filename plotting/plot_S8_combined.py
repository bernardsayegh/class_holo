import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
from getdist import loadMCSamples

# ============================================================
# CONFIGURATION
# ============================================================
CHAIN_ROOT = "/Users/Bernard_1/holographic_mcmc/chains_schw"
IGNORE_ROWS = 0.30

# Panel (a): bar chart models
MODELS = [
    (r'$\Lambda$CDM', "lcdm_shoes"),
    (r'Holographic'+"\n"+r'($\beta=1/12$)', "modelA_fixed"),
]

OBS = {
    "DES Y3": {"S8": 0.776, "S8_err": 0.017, "sigma8": 0.760, "sigma8_err": 0.020, "band": True},
    "KiDS-1000": {"S8": 0.759, "S8_err": 0.024, "sigma8": None, "sigma8_err": None, "band": False},
}

CLASS_THETAS = {
    r'$\Lambda$CDM': {"sigma8": 0.80761, "S8": 0.82841},
    r'Holographic'+"\n"+r'($\beta=1/12$)': {"sigma8": 0.74255, "S8": 0.77793},
}

# Panel (b): f_clust curve
# Endpoints from your MCMC:
#   f=0: S8 ~ 0.769 (Model A fixed)
#   f=1: S8 = 0.804 (Model A f=1)
FCLUST_GRID = np.linspace(0, 1, 11)
S8_THEORY_MED = np.linspace(0.769, 0.804, 11)  # Linear interpolation
S8_THEORY_LO = S8_THEORY_MED - 0.006
S8_THEORY_HI = S8_THEORY_MED + 0.006

# From Model A (f free): f_clust = 0.223 ± 0.154, S8 = 0.775 ± 0.008
MCMC_FREE_F = 0.223
MCMC_FREE_F_ERR = 0.154
MCMC_FREE_S8 = 0.775
MCMC_FREE_S8_ERR = 0.008

OUT_PNG = "S8_combined.png"
OUT_PDF = "S8_combined.pdf"

# ============================================================
# HELPERS
# ============================================================
def find_chain_base(chain_root, model_prefix):
    base = os.path.join(chain_root, model_prefix)
    if os.path.exists(base + ".1.txt"):
        return base
    matches = sorted(glob.glob(os.path.join(chain_root, model_prefix + "*.1.txt")))
    if matches:
        return matches[0].replace(".1.txt", "")
    raise FileNotFoundError(f"No chain files for '{model_prefix}' in '{chain_root}'")

def ensure_derived_S8(samples):
    names = [p.name for p in samples.paramNames.names]
    if "S8" in names:
        return
    if "sigma8" not in names or "Omega_m" not in names:
        raise ValueError("Cannot derive S8")
    i_sig8 = names.index("sigma8")
    i_Om = names.index("Omega_m")
    S8_vals = samples.samples[:, i_sig8] * np.sqrt(samples.samples[:, i_Om] / 0.3)
    samples.addDerived(S8_vals, name="S8", label="S_8")

def get_stats(samples, param):
    names = [p.name for p in samples.paramNames.names]
    i = names.index(param)
    m = samples.getMeans()
    cov = samples.getCov()
    return float(m[i]), float(np.sqrt(cov[i, i]))

# ============================================================
# LOAD MCMC
# ============================================================
results = {}
for label, prefix in MODELS:
    try:
        base = find_chain_base(CHAIN_ROOT, prefix)
        s = loadMCSamples(base, settings={"ignore_rows": IGNORE_ROWS})
        ensure_derived_S8(s)
        sig8, sig8_err = get_stats(s, "sigma8")
        S8, S8_err = get_stats(s, "S8")
        results[label] = dict(sigma8=sig8, sigma8_err=sig8_err, S8=S8, S8_err=S8_err)
        print(f"Loaded {label}: σ8={sig8:.4f}±{sig8_err:.4f}, S8={S8:.4f}±{S8_err:.4f}")
    except Exception as e:
        print(f"Warning: {label}: {e}, using CLASS values")
        cv = CLASS_THETAS.get(label, {})
        results[label] = dict(sigma8=cv.get("sigma8", 0.8), sigma8_err=0.006,
                              S8=cv.get("S8", 0.8), S8_err=0.005)

# Build arrays
labels_a = [m[0] for m in MODELS] + list(OBS.keys())
sig8_vals = [results[l]["sigma8"] for l in [m[0] for m in MODELS]]
sig8_errs = [results[l]["sigma8_err"] for l in [m[0] for m in MODELS]]
S8_vals = [results[l]["S8"] for l in [m[0] for m in MODELS]]
S8_errs = [results[l]["S8_err"] for l in [m[0] for m in MODELS]]

for k in OBS.keys():
    sig8_vals.append(OBS[k]["sigma8"] if OBS[k]["sigma8"] else np.nan)
    sig8_errs.append(OBS[k]["sigma8_err"] if OBS[k]["sigma8_err"] else 0.0)
    S8_vals.append(OBS[k]["S8"])
    S8_errs.append(OBS[k]["S8_err"])

# ============================================================
# CREATE FIGURE
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# -------------------- Panel (a): Bar chart --------------------
x = np.arange(len(labels_a))
w = 0.32

# DES band
for name, d in OBS.items():
    if d.get("band"):
        ax1.axhspan(d["S8"] - d["S8_err"], d["S8"] + d["S8_err"],
                    alpha=0.18, color="green", label=f"{name} band")

bars_sig8 = ax1.bar(x - w/2, sig8_vals, w, yerr=sig8_errs, capsize=5,
                     color='steelblue', label=r"$\sigma_8$")
bars_S8 = ax1.bar(x + w/2, S8_vals, w, yerr=S8_errs, capsize=5,
                   color='coral', label=r"$S_8$")

# Hide NaN bars
for i, v in enumerate(sig8_vals):
    if np.isnan(v):
        bars_sig8[i].set_alpha(0.0)

# Annotate
for b, v, e in zip(bars_sig8, sig8_vals, sig8_errs):
    if not np.isnan(v):
        ax1.text(b.get_x() + b.get_width()/2, v + e + 0.004, f"{v:.3f}",
                 ha="center", va="bottom", fontsize=10)
for b, v, e in zip(bars_S8, S8_vals, S8_errs):
    ax1.text(b.get_x() + b.get_width()/2, v + e + 0.004, f"{v:.3f}",
             ha="center", va="bottom", fontsize=10)

# CLASS dots
for i, lab in enumerate([m[0] for m in MODELS]):
    cv = CLASS_THETAS.get(lab)
    if cv:
        ax1.plot(i - w/2, cv["sigma8"], 'ko', ms=8, zorder=6)
        ax1.plot(i + w/2, cv["S8"], 'ko', ms=8, zorder=6)
ax1.plot([], [], 'ko', ms=8, label=r"CLASS ($\theta_s$-matched)")

ax1.set_title(r"(a) Resolution of the $S_8$ Tension", fontsize=14)
ax1.set_ylabel("Clustering Amplitude", fontsize=13)
ax1.set_xticks(x)
ax1.set_xticklabels(labels_a, fontsize=11)
ax1.set_ylim(0.70, 0.88)
ax1.legend(loc="upper right", fontsize=10, framealpha=0.95)

# -------------------- Panel (b): S8 vs f_clust --------------------
ax2.axhspan(0.776 - 0.017, 0.776 + 0.017, alpha=0.2, color="green")
ax2.axhline(0.776, color="green", ls="--", lw=2, label=r"DES Y3 $S_8 = 0.776 \pm 0.017$")
ax2.axhline(0.834, color="gray", ls="--", lw=1.5, label=r"Planck $\Lambda$CDM $S_8 = 0.834$")
ax2.axhline(0.805, color="red", ls=":", lw=1.5, label=r"MCMC $\Lambda$CDM $S_8 = 0.805$")

# Theory band (β = 1/12, varying f_clust)
ax2.fill_between(FCLUST_GRID, S8_THEORY_LO, S8_THEORY_HI, alpha=0.3, color="blue")
ax2.plot(FCLUST_GRID, S8_THEORY_MED, lw=2.5, color="blue", label=r"Holographic $\beta = 1/12$")

# MCMC free point (β, f_clust both free)
ax2.errorbar(MCMC_FREE_F, MCMC_FREE_S8,
             xerr=MCMC_FREE_F_ERR, yerr=MCMC_FREE_S8_ERR,
             fmt="s", ms=11, capsize=5, lw=2, color="red",
             label=r"MCMC ($\beta, f_{\rm clust}$ free)")

ax2.set_title(r"(b) $S_8$ vs Clustering Fraction", fontsize=14)
ax2.set_xlabel(r"$f_{\rm clust}$", fontsize=13)
ax2.set_ylabel(r"$S_8 = \sigma_8\sqrt{\Omega_m/0.3}$", fontsize=13)
ax2.set_xlim(-0.05, 1.05)
ax2.set_ylim(0.74, 0.86)
ax2.legend(loc="upper left", fontsize=10)
ax2.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=200, bbox_inches='tight')
plt.savefig(OUT_PDF, bbox_inches='tight')
print(f"\nSaved: {OUT_PNG} and {OUT_PDF}")
