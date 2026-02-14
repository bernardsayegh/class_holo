import sys
sys.path.insert(0, "python")
from classy import Class
import numpy as np

BASE = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    "h": 0.6766,
    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
    "tau_reio": 0.054,
    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "Omega_k": 0.0,
}

configs = {
    "ΛCDM (no mapping)": {},
    "acdm (Amap=2, β=0)": {
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 0.0,
        "super_schw_Amap": 2.0,
        "super_schw_deltaS": 0.03,
        "super_schw_gamma": 2.0,
    },
    "Model A (β=1/12)": {
        "interaction_beta": 1.0/12.0,
        "f_clust": 0.0,
    },
    "Model B (Amap=2, β=1/12)": {
        "interaction_beta": 1.0/12.0,
        "f_clust": 0.0,
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 0.0,
        "super_schw_Amap": 2.0,
        "super_schw_deltaS": 0.03,
        "super_schw_gamma": 2.0,
    },
}

header = f"{'Model':<28} {'H0':>7} {'H0_loc':>8} {'sig8':>8} {'S8':>8} {'S8_loc':>8} {'Om':>8} {'Om_loc':>8} {'X0':>8}"
print(header)
print("=" * len(header))

for name, extra in configs.items():
    p = dict(BASE)
    p.update(extra)

    c = Class()
    c.set(p)
    c.compute()

    H0 = 100.0 * c.h()
    Om = float(c.Omega_m())
    sig8 = float(c.sigma8())
    S8 = sig8 * np.sqrt(Om / 0.3)

    try:
        d = c.get_current_derived_parameters(["H0_local", "X0_schw"])
        H0_loc = d.get("H0_local", H0)
        X0 = d.get("X0_schw", 0.0)
    except:
        H0_loc = H0
        X0 = 0.0

    # Local-frame quantities
    omega_m = Om * (H0 / 100.0) ** 2   # physical omega_m is fixed
    Om_loc = omega_m / (H0_loc / 100.0) ** 2
    S8_loc = sig8 * np.sqrt(Om_loc / 0.3)

    print(f"{name:<28} {H0:>7.2f} {H0_loc:>8.2f} {sig8:>8.5f} {S8:>8.5f} {S8_loc:>8.5f} {Om:>8.5f} {Om_loc:>8.5f} {X0:>8.5f}")

    c.struct_cleanup()
    c.empty()

print()
print("Key points:")
print("  - acdm: Amap=2 but β=0 → S never exceeds 1 → X0=0 → H0_loc=H0 → mapping does nothing")
print("  - Model A: β=1/12 but no mapping → σ8 suppressed, H0_loc≈H0")
print("  - Model B: β=1/12 + Amap=2 → σ8 suppressed AND H0_loc boosted")
print("  - S8_loc < S8 for Model B because Om_loc < Om (same ωm, larger H0_loc)")
