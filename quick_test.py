import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np

print("=" * 55)
print("HOLOGRAPHIC DARK DRAG VALIDATION")
print("=" * 55)

common = {
    'output': 'mPk',
    'P_k_max_1/Mpc': 10.0,
    'h': 0.6766,
    'omega_b': 0.02242,
    'omega_cdm': 0.11933,
    'n_s': 0.9665,
    'A_s': 2.1e-9,
}

print(f"\n{'β':>8} | {'σ₈':>10} | {'Ω_m':>10} | {'S₈':>10}")
print("-" * 50)

for beta in [0.0, -0.25, -0.50, -0.75]:
    cosmo = Class()
    params = common.copy()
    params['interaction_beta'] = beta
    cosmo.set(params)
    try:
        cosmo.compute()
        sig8 = cosmo.sigma8()
        Om = cosmo.Omega_m()
        S8 = sig8 * np.sqrt(Om / 0.3)
        marker = " ✓" if abs(S8 - 0.770) < 0.02 else ""
        print(f"{beta:>8.2f} | {sig8:>10.4f} | {Om:>10.4f} | {S8:>10.4f}{marker}")
    except Exception as e:
        print(f"{beta:>8.2f} | Error: {e}")
    cosmo.struct_cleanup()
    cosmo.empty()

print("-" * 50)
print("Target: S₈ = 0.770 (weak lensing)")
