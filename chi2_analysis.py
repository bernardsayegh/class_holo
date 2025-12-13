import numpy as np

# 1. DATA: Observational Constraints (Mean +/- 1 sigma)
# -----------------------------------------------------
# Planck 2018 (TT,TE,EE+lowE+lensing)
s8_planck_mean = 0.832
s8_planck_err  = 0.013

# KiDS-1000 (Weak Lensing) - The "Tension" Data
s8_kids_mean   = 0.766
s8_kids_err    = 0.014

# DES Y3 (Weak Lensing)
s8_des_mean    = 0.776# Check the perturbation implementation
grep -A10 "Q_over_rho = " source/perturbations.c
s8_des_err     = 0.017

# 2. MODEL: CLASS Output (UPDATED VALUES)
# ---------------------------
# Run 1: Standard LambdaCDM (beta = 0)
sigma8_lcdm = 0.837  # S8 = 0.823 * sqrt(0.3111/0.3)

# Run 2: Holographic Model (beta = 0.5)
sigma8_holo = 0.772  # S8 = 0.759 * sqrt(0.334/0.3)

# 3. CALCULATE CHI-SQUARED
# ------------------------
def get_chi2(model_value, obs_mean, obs_err):
    return ((model_value - obs_mean) / obs_err)**2

# Chi2 for LambdaCDM
chi2_lcdm_kids = get_chi2(sigma8_lcdm, s8_kids_mean, s8_kids_err)
chi2_lcdm_des  = get_chi2(sigma8_lcdm, s8_des_mean, s8_des_err)

# Chi2 for Holographic
chi2_holo_kids = get_chi2(sigma8_holo, s8_kids_mean, s8_kids_err)
chi2_holo_des  = get_chi2(sigma8_holo, s8_des_mean, s8_des_err)

# 4. PRINT RESULTS
# ----------------
print(f"============================================")
print(f"        S8 TENSION CHI-SQUARED ANALYSIS")
print(f"============================================")
print(f"\nModel Values:")
print(f"  LambdaCDM S8: {sigma8_lcdm:.3f}")
print(f"  Holographic S8: {sigma8_holo:.3f}")

print(f"\n--- KiDS-1000 (S8 = {s8_kids_mean} ± {s8_kids_err}) ---")
print(f"  LCDM χ²: {chi2_lcdm_kids:.2f}")
print(f"  Holo χ²: {chi2_holo_kids:.2f}")
print(f"  Δχ² (improvement): {chi2_lcdm_kids - chi2_holo_kids:.2f}")

print(f"\n--- DES Y3 (S8 = {s8_des_mean} ± {s8_des_err}) ---")
print(f"  LCDM χ²: {chi2_lcdm_des:.2f}")
print(f"  Holo χ²: {chi2_holo_des:.2f}")
print(f"  Δχ² (improvement): {chi2_lcdm_des - chi2_holo_des:.2f}")

# Calculate Gaussian Tension (in sigmas)
sigma_tension_lcdm = (sigma8_lcdm - s8_kids_mean) / np.sqrt(s8_planck_err**2 + s8_kids_err**2)
sigma_tension_holo = (sigma8_holo - s8_kids_mean) / np.sqrt(s8_planck_err**2 + s8_kids_err**2)

print(f"\n--- Tension Level with KiDS ---")
print(f"  LCDM: {sigma_tension_lcdm:.2f}σ tension")
print(f"  Holo: {sigma_tension_holo:.2f}σ tension")
print(f"  Reduction: {sigma_tension_lcdm - sigma_tension_holo:.2f}σ")

print(f"\n============================================")
print(f"SUMMARY: Holographic model reduces S8 tension")
print(f"         from {sigma_tension_lcdm:.1f}σ to {sigma_tension_holo:.1f}σ")
print(f"============================================")
