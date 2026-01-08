import numpy as np
import sys
sys.path.insert(0, "python")
from classy import Class

def run_model(name, params):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    
    sigma8 = cosmo.sigma8()
    Omega_m = cosmo.Omega_m()
    S8 = sigma8 * np.sqrt(Omega_m / 0.3)
    h = cosmo.h()
    H0 = 100 * h
    
    print(f"{name}:")
    print(f"  sigma8  = {sigma8:.4f}")
    print(f"  Omega_m = {Omega_m:.4f}")
    print(f"  S8      = {S8:.4f}")
    print(f"  H0      = {H0:.2f}")
    
    cosmo.struct_cleanup()
    cosmo.empty()
    return sigma8, Omega_m, S8, H0

base = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    "z_max_pk": 1.0,
    "h": 0.6766,
    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
}

# LCDM
lcdm = {**base, "Omega_Lambda": 0.6847}

# Holographic beta = 0.5 (original)
holo_05 = {
    **base,
    "Omega_Lambda": 0.0,
    "w0_fld": -1.0,
    "wa_fld": 0.0,
    "interaction_beta": 0.5,
    "interaction_area_dilution": 1,
    "interaction_use_ah_filter": 1,
    "f_clust": 0.0,
}

# Holographic beta = 0.25 (half strength)
holo_025 = {
    **base,
    "Omega_Lambda": 0.0,
    "w0_fld": -1.0,
    "wa_fld": 0.0,
    "interaction_beta": 0.25,
    "interaction_area_dilution": 1,
    "interaction_use_ah_filter": 1,
    "f_clust": 0.0,
}

# Holographic beta = 0.125 (effective coupling directly)
holo_0125 = {
    **base,
    "Omega_Lambda": 0.0,
    "w0_fld": -1.0,
    "wa_fld": 0.0,
    "interaction_beta": 0.125,
    "interaction_area_dilution": 1,
    "interaction_use_ah_filter": 1,
    "f_clust": 0.0,
}

print("="*50)
run_model("LCDM", lcdm)
print()
run_model("Holographic beta=0.50 (beta_eff=0.125)", holo_05)
print()
run_model("Holographic beta=0.25 (beta_eff=0.0625)", holo_025)
print()
run_model("Holographic beta=0.125 (beta_eff=0.03125)", holo_0125)
print("="*50)

print("\nComparison to observations:")
print("  KiDS-1000: S8 = 0.766 +/- 0.014")
print("  DES Y3:    S8 = 0.776 +/- 0.017")
print("  Planck:    S8 = 0.834 +/- 0.013")
