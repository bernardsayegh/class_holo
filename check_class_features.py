#!/usr/bin/env python3
import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np

C_KMS = 299792.458
THETA_S_TARGET = 1.040423

BASE = {
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

SCR_MAP = {
    "super_schwarzschild_correction": "yes",
    "super_schw_amp": 0.0,
    "super_schw_Amap": 2.0,
    "super_schw_deltaS": 0.03,
    "super_schw_gamma": 2.0,
    "super_schw_no_mapping": 0,
}

def run(label, extra):
    params = dict(BASE)
    params.update(extra)
    def res(h):
        p = dict(params); p["h"] = float(h)
        c = Class(); c.set(p); c.compute()
        t = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
        c.struct_cleanup(); c.empty()
        return t - THETA_S_TARGET
    h = brentq(res, 0.55, 0.80, xtol=1e-8)
    params["h"] = float(h)
    c = Class(); c.set(params); c.compute()
    H0 = c.Hubble(0) * C_KMS
    s8 = c.sigma8()
    Om = c.Omega_m()
    S8 = s8 * np.sqrt(Om/0.3)
    H0_loc = c.get_current_derived_parameters(["H0_local"])["H0_local"]
    print(f"  {label}:")
    print(f"    H0={H0:.2f}  H0_local={H0_loc:.2f}  sigma8={s8:.4f}  S8={S8:.4f}  Om={Om:.4f}")
    c.struct_cleanup(); c.empty()

print("=" * 60)
run("LCDM", {})
run("LCDM + mapping (beta=0, Amap=2)", SCR_MAP)
run("Model B (beta=1/12, Amap=2)", dict(SCR_MAP, interaction_beta=1./12., f_clust=0.0))
print("=" * 60)
