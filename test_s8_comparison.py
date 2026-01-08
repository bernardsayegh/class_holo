import numpy as np
import sys
sys.path.insert(0, "python")
from classy import Class

def get_s8(params):
    cosmo = Class()
    full_params = {
        "output": "mPk",
        "P_k_max_1/Mpc": 10.0,
        "z_max_pk": 1.0,
    }
    full_params.update(params)
    
    # Convert ln10^10A_s to A_s if needed
    if 'ln10^{10}A_s' in full_params:
        ln10_10_As = full_params.pop('ln10^{10}A_s')
        full_params['A_s'] = np.exp(ln10_10_As) * 1e-10
    
    cosmo.set(full_params)
    cosmo.compute()
    sigma8 = cosmo.sigma8()
    Omega_m = cosmo.Omega_m()
    s8 = sigma8 * np.sqrt(Omega_m / 0.3)
    
    print(f"    sigma8={sigma8:.4f}, Omega_m={Omega_m:.4f}, H0={cosmo.h()*100:.2f}")
    
    cosmo.struct_cleanup()
    cosmo.empty()
    return s8

# 1. Standard Planck LCDM Baseline
params_lcdm = {
    'omega_b': 0.02237,
    'omega_cdm': 0.1200, 
    'h': 0.6736,
    'ln10^{10}A_s': 3.044,
    'n_s': 0.9649,
    'tau_reio': 0.0544,
    'interaction_beta': 0.0,
    'Omega_Lambda': 0.6847,
}

# 2. Your Python Test (Beta=0.22, but old parameters)
params_mixed = {
    'omega_b': 0.02237,
    'omega_cdm': 0.1200, 
    'h': 0.6736,
    'ln10^{10}A_s': 3.044,
    'n_s': 0.9649,
    'tau_reio': 0.0544,
    'interaction_beta': 0.22,
    'Omega_Lambda': 0.0,
    'w0_fld': -1.0,
    'wa_fld': 0.0,
    'interaction_area_dilution': 1,
    'interaction_use_ah_filter': 1,
    'f_clust': 0.0,
}

# 3. The ACTUAL MCMC Best-Fit (From your logs)
params_mcmc = {
    'omega_b': 0.02244,
    'omega_cdm': 0.1175,
    'h': 0.6809,
    'ln10^{10}A_s': 3.013,
    'n_s': 0.968,
    'tau_reio': 0.043,
    'interaction_beta': 0.218,
    'Omega_Lambda': 0.0,
    'w0_fld': -1.0,
    'wa_fld': 0.0,
    'interaction_area_dilution': 1,
    'interaction_use_ah_filter': 1,
    'f_clust': 0.0,
}

print("="*60)
print("S8 COMPARISON")
print("="*60)

print("\n1. LCDM Baseline (Planck 2018):")
s8_lcdm = get_s8(params_lcdm)
print(f"   S8 = {s8_lcdm:.4f}")

print("\n2. Beta=0.22 with LCDM parameters (mixed):")
s8_mixed = get_s8(params_mixed)
print(f"   S8 = {s8_mixed:.4f}")

print("\n3. Full MCMC Best-Fit (beta=0.218):")
s8_mcmc = get_s8(params_mcmc)
print(f"   S8 = {s8_mcmc:.4f}")

print("\n" + "="*60)
print("OBSERVATIONS:")
print("  KiDS-1000: S8 = 0.766 +/- 0.014")
print("  DES Y3:    S8 = 0.776 +/- 0.017")
print("  Planck:    S8 = 0.834 +/- 0.013")
print("="*60)
