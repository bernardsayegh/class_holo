import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from getdist import plots, loadMCSamples
import os

# Target values with errors
H0_SHOES = 73.04
H0_SHOES_ERR = 1.04
S8_DES = 0.776
S8_DES_ERR = 0.017
OMEGA_M_PLANCK = 0.315
OMEGA_M_PLANCK_ERR = 0.007

SIGMA8_WL = S8_DES / np.sqrt(OMEGA_M_PLANCK / 0.3)
SIGMA8_WL_ERR = S8_DES_ERR / np.sqrt(OMEGA_M_PLANCK / 0.3)

print(f"Derived σ8 from S8: {SIGMA8_WL:.3f} ± {SIGMA8_WL_ERR:.3f}")

# Load chains
samples_list = []
configs = [
    ("chains_schw/lcdm_shoes", "ΛCDM", "gray", "H0"),
    ("chains_schw/modelA_fixed", "Model A (β=1/12)", "blue", "H0"),
    ("chains_schw/modelB_Amap2", "Model B (Amap=2)", "red", "H0_local"),
]

for path, label, color, h0_param in configs:
    if os.path.exists(f"{path}.1.txt"):
        s = loadMCSamples(path, settings={"ignore_rows": 0.3})
        s.label = label
        s.color = color
        param_names = [p.name for p in s.paramNames.names]
        
        if "S8" not in param_names:
            s8_idx = param_names.index("sigma8")
            Om_idx = param_names.index("Omega_m")
            S8_vals = s.samples[:, s8_idx] * np.sqrt(s.samples[:, Om_idx] / 0.3)
            s.addDerived(S8_vals, name="S8", label="S_8")
        
        if h0_param == "H0_local" and "H0_local" in param_names:
            h0_idx = param_names.index("H0_local")
        else:
            h0_idx = param_names.index("H0")
        H0_eff_vals = s.samples[:, h0_idx]
        s.addDerived(H0_eff_vals, name="H0_eff", label="H_0^{\\rm eff}")
        
        samples_list.append(s)
        print(f"Loaded {label}: N={s.numrows}, using {h0_param}")
    else:
        print(f"Missing: {path}")

params = ["H0_eff", "sigma8", "Omega_m", "S8"]

g = plots.get_subplot_plotter(width_inch=12)
g.settings.alpha_filled_add = 0.7
g.settings.axes_fontsize = 13
g.settings.lab_fontsize = 15
g.settings.legend_fontsize = 11
g.settings.tight_layout = True
g.settings.subplot_size_ratio = 0.85

g.triangle_plot(
    samples_list,
    params,
    filled=True,
    contour_colors=["gray", "blue", "red"],
    legend_loc=None,
)

targets = {
    "H_0": (H0_SHOES, H0_SHOES_ERR, "green"),
    "sigma": (SIGMA8_WL, SIGMA8_WL_ERR, "purple"),
    "Omega": (OMEGA_M_PLANCK, OMEGA_M_PLANCK_ERR, "cyan"),
    "S_8": (S8_DES, S8_DES_ERR, "orange"),
}

for ax in g.subplots.flatten():
    if ax is None:
        continue
    xlabel = ax.get_xlabel()
    ylabel = ax.get_ylabel()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    for key, (val, err, color) in targets.items():
        if key in xlabel:
            ax.axvline(val, color=color, ls="-", lw=1.5, alpha=0.8)
            ax.axvspan(val - err, val + err, color=color, alpha=0.15, zorder=0)
        if key in ylabel:
            ax.axhline(val, color=color, ls="-", lw=1.5, alpha=0.8)
            ax.axhspan(val - err, val + err, color=color, alpha=0.15, zorder=0)
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

# Remove any auto-generated legends
for ax in g.subplots.flatten():
    if ax is None:
        continue
    lg = ax.get_legend()
    if lg is not None:
        lg.remove()
for lg in list(g.fig.legends):
    lg.remove()

# Custom legends
model_handles = [
    Patch(facecolor="gray", alpha=0.5, edgecolor="gray", label="ΛCDM"),
    Patch(facecolor="blue", alpha=0.5, edgecolor="blue", label="Model A (β=1/12)"),
    Patch(facecolor="red",  alpha=0.5, edgecolor="red",  label="Model B ($A_{\\rm map}$=2)"),
]

obs_handles = [
    Patch(facecolor="green",  alpha=0.3, edgecolor="green",  label="SH0ES ($H_0$)"),
    Patch(facecolor="purple", alpha=0.3, edgecolor="purple", label="DES-Y3 ($\\sigma_8$)"),
    Patch(facecolor="cyan",   alpha=0.3, edgecolor="cyan",   label="Planck ($\\Omega_m$)"),
    Patch(facecolor="orange", alpha=0.3, edgecolor="orange", label="DES-Y3 ($S_8$)"),
]

leg_models = g.fig.legend(
    handles=model_handles,
    loc="upper right",
    bbox_to_anchor=(0.97, 0.97),
    fontsize=10,
    title="Models",
    title_fontsize=11,
    frameon=True,
    fancybox=True,
)

leg_constraints = g.fig.legend(
    handles=obs_handles,
    loc="upper right",
    bbox_to_anchor=(0.97, 0.80),
    fontsize=10,
    title="Constraints",
    title_fontsize=11,
    frameon=True,
    fancybox=True,
)

g.fig.add_artist(leg_models)

g.fig.suptitle("Holographic Cosmology: Resolving H₀ and S₈ Tensions", fontsize=15, y=1.01)

plt.savefig("triangle_tensions.pdf", bbox_inches="tight", dpi=300)
plt.savefig("triangle_tensions.png", bbox_inches="tight", dpi=300)
print("\nSaved: triangle_tensions.pdf and triangle_tensions.png")
plt.show()
