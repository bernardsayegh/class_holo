#!/usr/bin/env python3
"""Corner plot H0_eff/sigma8/Omega_m/S8 with labeled constraint bands
(fig:triangle). Run from anywhere."""
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from getdist import plots, loadMCSamples

CHAINS_DIR = os.path.expanduser("~/class_holo_test/chains_test")
H0_SHOES, H0_SHOES_ERR = 73.49, 0.93          # SH0ES R25
S8_DES, S8_DES_ERR = 0.776, 0.017             # DES-Y3 WL
S8_GROWTH, S8_GROWTH_ERR = 0.747, 0.029       # fsigma8 growth (Nguyen+ 2023, PRL 131 111001)
OM_PLANCK, OM_PLANCK_ERR = 0.315, 0.007
OM_PANTHEON, OM_PANTHEON_ERR = 0.334, 0.018   # Pantheon+ shape

def find_root(candidates):
    for c in candidates:
        if os.path.exists(os.path.join(CHAINS_DIR, c + ".1.txt")):
            return os.path.join(CHAINS_DIR, c)
    hits = sorted(glob.glob(os.path.join(CHAINS_DIR, candidates[0].split("_v5")[0] + "*v5*.1.txt")))
    return hits[0][:-6] if hits else None

configs = [
    (["lcdm_v5_acc", "lcdm_v5_plancknu", "lcdm_v5"], "LCDM", "gray", "H0"),
    (["modelA_v5_acc", "modelA_v5"], "Model A", "blue", "H0"),
    (["modelB_A2_v5_acc", "modelB_A2_v5"], "Model B", "red", "H0_local"),
]
samples_list = []
for candidates, label, color, h0_param in configs:
    root = find_root(candidates)
    if root is None:
        print(f"Missing v5 chain for {candidates[0]}"); continue
    s = loadMCSamples(root, settings={"ignore_rows": 0.3})
    s.label, s.color = label, color
    pn = [p.name for p in s.paramNames.names]
    if "S8" not in pn:
        S8v = s.samples[:, pn.index("sigma8")] * np.sqrt(s.samples[:, pn.index("Omega_m")] / 0.3)
        s.addDerived(S8v, name="S8", label="S_8")
    h0_idx = pn.index("H0_local") if (h0_param == "H0_local" and "H0_local" in pn) else pn.index("H0")
    s.addDerived(s.samples[:, h0_idx], name="H0_eff", label="H_0^{\\rm eff}")
    w = s.weights
    sig8 = s.samples[:, pn.index("sigma8")]
    Om = s.samples[:, pn.index("Omega_m")]
    S8 = sig8 * np.sqrt(Om / 0.3)
    H0 = s.samples[:, h0_idx]
    print(f"Loaded {label}: {os.path.basename(root)}  {h0_param}={np.average(H0,weights=w):.2f}, "
          f"S8={np.average(S8,weights=w):.3f}, sig8={np.average(sig8,weights=w):.3f}, "
          f"Om={np.average(Om,weights=w):.3f}")
    samples_list.append(s)

params = ["H0_eff", "sigma8", "Omega_m", "S8"]
g = plots.get_subplot_plotter(width_inch=12)
g.settings.alpha_filled_add = 0.7
g.settings.axes_fontsize = 13
g.settings.lab_fontsize = 15
g.settings.legend_fontsize = 11
g.settings.tight_layout = True
g.settings.subplot_size_ratio = 0.85
g.triangle_plot(samples_list, params, filled=True,
                contour_colors=["gray", "blue", "red"], legend_loc=None)

targets = {
    "H_0":   [(H0_SHOES, H0_SHOES_ERR, "green")],
    "Omega": [(OM_PLANCK, OM_PLANCK_ERR, "cyan"),
              (OM_PANTHEON, OM_PANTHEON_ERR, "olive")],
    "S_8":   [(S8_DES, S8_DES_ERR, "orange"),
              (S8_GROWTH, S8_GROWTH_ERR, "purple")],
}
for ax in g.subplots.flatten():
    if ax is None: continue
    xlabel, ylabel = ax.get_xlabel(), ax.get_ylabel()
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    for key, bands in targets.items():
        for val, err, color in bands:
            if key in xlabel:
                ax.axvline(val, color=color, ls="-", lw=1.5, alpha=0.8)
                ax.axvspan(val - err, val + err, color=color, alpha=0.15, zorder=0)
            if key in ylabel:
                ax.axhline(val, color=color, ls="-", lw=1.5, alpha=0.8)
                ax.axhspan(val - err, val + err, color=color, alpha=0.15, zorder=0)
    ax.set_xlim(xlim); ax.set_ylim(ylim)

for ax in g.subplots.flatten():
    if ax is None: continue
    lg = ax.get_legend()
    if lg is not None: lg.remove()
for lg in list(g.fig.legends): lg.remove()

model_handles = [
    Patch(facecolor="gray", alpha=0.5, edgecolor="gray", label="LCDM"),
    Patch(facecolor="blue", alpha=0.5, edgecolor="blue", label=r"Model A ($\beta=1/12$)"),
    Patch(facecolor="red",  alpha=0.5, edgecolor="red",  label=r"Model B ($A_{\rm map}=2$)"),
]
obs_handles = [
    Patch(facecolor="green",  alpha=0.3, edgecolor="green",  label=r"SH0ES R25 ($H_0$)"),
    Patch(facecolor="orange", alpha=0.3, edgecolor="orange", label=r"DES-Y3 WL ($S_8$)"),
    Patch(facecolor="purple", alpha=0.3, edgecolor="purple", label=r"Growth $f\sigma_8$ ($S_8$)"),
    Patch(facecolor="cyan",   alpha=0.3, edgecolor="cyan",   label=r"Planck ($\Omega_m$)"),
    Patch(facecolor="olive",  alpha=0.3, edgecolor="olive",  label=r"Pantheon+ shape ($\Omega_m$)"),
]
leg_models = g.fig.legend(handles=model_handles, loc="upper right",
                          bbox_to_anchor=(0.97, 0.97), fontsize=10,
                          title="Models", title_fontsize=11, frameon=True, fancybox=True)
g.fig.legend(handles=obs_handles, loc="upper right",
             bbox_to_anchor=(0.97, 0.80), fontsize=10,
             title="Constraints", title_fontsize=11, frameon=True, fancybox=True)
g.fig.add_artist(leg_models)
g.fig.suptitle("Holographic Cosmology: Resolving H0 and S8 Tensions", fontsize=15, y=1.01)
plt.savefig("triangle_S8_Om_H0.pdf", bbox_inches="tight", dpi=300)
plt.savefig("triangle_S8_Om_H0.png", bbox_inches="tight", dpi=200)
print(f"\nSaved: {os.path.abspath('triangle_S8_Om_H0.pdf')}  [matches includegraphics name]")
