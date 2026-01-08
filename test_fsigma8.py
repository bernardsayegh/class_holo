import numpy as np
import sys
sys.path.insert(0, "python")
from classy import Class

z_obs = np.array([0.38, 0.51, 0.61])
fs8_obs = np.array([0.497, 0.458, 0.436])
sig_obs = np.array([0.045, 0.038, 0.034])

def get_fsigma8(cosmo, z):
    dz = 1e-3
    s0 = float(cosmo.sigma(8.0 / cosmo.h(), z))
    sp = float(cosmo.sigma(8.0 / cosmo.h(), z + dz))
    sm = float(cosmo.sigma(8.0 / cosmo.h(), max(z - dz, 0.0)))
    dlns_dz = (np.log(sp) - np.log(sm)) / (2.0 * dz)
    f = -(1.0 + z) * dlns_dz
    return f * s0

def chi2_diag(model):
    return float(np.sum(((model - fs8_obs) / sig_obs) ** 2))

def run_model(name, params):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    fs8 = np.array([get_fsigma8(cosmo, z) for z in z_obs])
    c2 = chi2_diag(fs8)
    print(f"{name}")
    for i, z in enumerate(z_obs):
        print(f"  z={z:.2f}  fs8={fs8[i]:.4f}  (obs: {fs8_obs[i]:.3f})")
    print(f"  chi2 = {c2:.2f}")
    cosmo.struct_cleanup()
    cosmo.empty()
    return fs8, c2

lcdm_params = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    "z_max_pk": 1.0,
    "h": 0.6766,
    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
    "Omega_Lambda": 0.6847,
}

holo_params = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    "z_max_pk": 1.0,
    "h": 0.6766,
    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
    "Omega_Lambda": 0.0,
    "w0_fld": -1.0,
    "wa_fld": 0.0,
    "interaction_beta": 0.5,
    "interaction_area_dilution": 1,
    "interaction_use_ah_filter": 1,
    "f_clust": 0.0,
}

fs8_l, c2_l = run_model("LCDM", lcdm_params)
print()
fs8_h, c2_h = run_model("Holographic", holo_params)
print()
print(f"Delta chi2 (Holo - LCDM) = {c2_h - c2_l:.2f}")
