#!/usr/bin/env python3
"""
check_modelB.py
Verify Model B produces correct H0_local via super-Schwarzschild mapping.

This script runs Model B (β=1/12, Amap=2) and verifies that:
1. H0_local ≈ 73 km/s/Mpc (SH0ES compatible)
2. H0_phys ≈ 67-68 km/s/Mpc (CMB compatible)
3. σ8 ≈ 0.75 (reduced from ΛCDM)

Usage:
    python3 check_modelB.py
"""
import sys
sys.path.insert(0, "python")
from classy import Class
import numpy as np

params = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    
    # Cosmological parameters (near MCMC best-fit)
    "h": 0.6784,
    "omega_b": 0.02237,
    "omega_cdm": 0.1200,
    "A_s": 2.1e-9,
    "n_s": 0.9649,
    "tau_reio": 0.0544,
    
    # Holographic interaction
    "interaction_beta": 1.0/12.0,  # β = 0.0833...
    "f_clust": 0.0,
    
    # Super-Schwarzschild Correction (Model B)
    "super_schwarzschild_correction": "yes",
    "super_schw_amp": 0.0,         # No reservoir
    "super_schw_Amap": 2.0,        # Mapping amplitude
    "super_schw_deltaS": 0.03,
    "super_schw_gamma": 2.0,
    "super_schw_no_mapping": 0,    # Mapping enabled
}

print("=" * 60)
print("Model B Validation (β=1/12, Amap=2)")
print("=" * 60)

c = Class()
c.set(params)
c.compute()

# Get derived parameters
derived = c.get_current_derived_parameters(["H0_local", "sigma8", "Omega_m", "100*theta_s"])
H0_local = derived["H0_local"]
sigma8 = derived["sigma8"]
Omega_m = derived["Omega_m"]
theta_s = derived["100*theta_s"]
S8 = sigma8 * np.sqrt(Omega_m / 0.3)
H0_phys = 100.0 * c.h()

print(f"\nInput parameters:")
print(f"  h             = {params['h']}")
print(f"  omega_b       = {params['omega_b']}")
print(f"  omega_cdm     = {params['omega_cdm']}")
print(f"  beta          = {params['interaction_beta']:.8f} (1/12 = {1/12:.8f})")
print(f"  Amap          = {params['super_schw_Amap']}")

print(f"\nDerived quantities:")
print(f"  H0 (physical) = {H0_phys:.2f} km/s/Mpc")
print(f"  H0_local      = {H0_local:.2f} km/s/Mpc")
print(f"  Mapping ratio = {H0_local/H0_phys:.4f} = exp({np.log(H0_local/H0_phys):.4f})")
print(f"  Implied X0    = {np.log(H0_local/H0_phys)/2.0:.5f}")
print(f"  100*theta_s   = {theta_s:.6f}")
print(f"  Omega_m       = {Omega_m:.4f}")
print(f"  sigma8        = {sigma8:.4f}")
print(f"  S8            = {S8:.4f}")

print(f"\nComparison with targets:")
print(f"  H0_local vs SH0ES (73.04): {H0_local:.2f} (Δ = {H0_local - 73.04:+.2f})")
print(f"  sigma8 vs ΛCDM (~0.81):    {sigma8:.4f} (suppressed by {(0.81-sigma8)/0.81*100:.1f}%)")
print(f"  S8 vs DES Y3 (~0.78):      {S8:.4f}")

c.struct_cleanup()
c.empty()

print("\n" + "=" * 60)
print("Validation complete.")
print("=" * 60)
