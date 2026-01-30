#!/usr/bin/env python3
import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np
C_KMS = 299792.458
THETA_S_TARGET = 1.040423
H_BRACKET = (0.55, 0.80)
BETA = 1.0/12.0
BASE_PARAMS = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
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
MODELS = {
    "LCDM": {},
    "Model A": {"interaction_beta": BETA, "f_clust": 0.0},
    "Model B": {"interaction_beta": BETA, "f_clust": 0.0,
                "super_schwarzschild_correction": "yes", "super_schw_amp": 0.0,
                "super_schw_Amap": 2.0, "super_schw_deltaS": 0.03,
                "super_schw_gamma": 2.0, "super_schw_no_mapping": 0},
    "Model C": {"interaction_beta": BETA, "f_clust": 0.0,
                "super_schwarzschild_correction": "yes", "super_schw_amp": 1.0,
                "super_schw_Amap": 2.0, "super_schw_deltaS": 0.03,
                "super_schw_gamma": 2.0, "super_schw_no_mapping": 0},
    "Model D": {"interaction_beta": BETA, "f_clust": 0.0,
                "super_schwarzschild_correction": "yes", "super_schw_amp": 1.0,
                "super_schw_Amap": 2.0, "super_schw_deltaS": 0.03,
                "super_schw_gamma": 2.0, "super_schw_no_mapping": 1},
}
def H0_phys_from_bg(c: Class) -> float:
    bg = c.get_background()
    z = bg["z"]
    i0 = int(np.argmin(np.abs(z)))
    return float(bg["H [1/Mpc]"][i0] * C_KMS)
def compute_stats(h, model_params):
    p = dict(BASE_PARAMS)
    p["h"] = float(h)
    p.update(model_params)
    c = Class()
    c.set(p)
    c.compute()
    theta = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
    H0_in = 100.0 * c.h()
    H0_phys = H0_phys_from_bg(c)
    Om = float(c.Omega_m())
    sig8 = float(c.sigma8())
    S8 = sig8 * np.sqrt(Om / 0.3)
    H0_loc = float(c.get_current_derived_parameters(["H0_local"])["H0_local"])
    X0 = 0.0
    try:
        X0 = float(c.get_current_derived_parameters(["X0_schw"])["X0_schw"])
    except Exception:
        pass
    c.struct_cleanup()
    c.empty()
    return dict(theta=float(theta), H0_in=H0_in, H0_phys=H0_phys, H0_loc=H0_loc,
                X0=X0, Omega_m=Om, sigma8=sig8, S8=S8)
def solve_h(model_params):
    def f(h):
        return compute_stats(h, model_params)["theta"] - THETA_S_TARGET
    return brentq(f, H_BRACKET[0], H_BRACKET[1], xtol=1e-8)
results = {}
for name, params in MODELS.items():
    print(f"Running {name}...")
    h = solve_h(params)
    results[name] = {"h": h, **compute_stats(h, params)}
print("\nSUMMARY (explicit H0_in vs H0_phys)")
print(f"{'Model':<10} {'h':>7} {'H0_in':>8} {'H0_phy':>8} {'H0_loc':>8} {'X0':>8} {'Om':>8} {'sig8':>8} {'S8':>8}")
for name, s in results.items():
    print(f"{name:<10} {s['h']:>7.5f} {s['H0_in']:>8.2f} {s['H0_phys']:>8.2f} {s['H0_loc']:>8.2f} "
          f"{s['X0']:>8.5f} {s['Omega_m']:>8.5f} {s['sigma8']:>8.5f} {s['S8']:>8.5f}")
print("\nTargets: SH0ES H0=73.04, DES S8=0.776")
