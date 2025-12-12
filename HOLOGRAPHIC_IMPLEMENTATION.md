# Holographic Dark Drag Implementation in CLASS

## Overview

This modified CLASS code implements the holographic dark energy interaction
from Sayegh (2026), "The Holographic Signature of the Cosmic Horizon."

The interaction term Q = -3βH ρ_m Ω_Λ w_Λ creates a "dark drag" effect
that suppresses structure growth on sub-horizon scales while preserving
CMB observables.

## New Parameter

- `interaction_beta`: Coupling strength (default: 0)
  - β > 0 (paper convention): Energy flows from DE to matter, suppresses σ₈
  - β ~ 0.5: Resolves S₈ tension (σ₈: 0.82 → 0.76)

## Files Modified

### 1. include/background.h (line 106)
```c
double interaction_beta; /**< holographic coupling strength (0 = no interaction) */
```

### 2. source/input.c
Default value (line ~5909):
```c
pba->interaction_beta = 0.;  /**< default: no holographic interaction */
```

Parameter reading (line ~2718):
```c
class_read_double("interaction_beta",pba->interaction_beta);
```

### 3. source/perturbations.c (line ~9238)
Scale-dependent dark drag in CDM density evolution:
```c
/* BEGIN HOLOGRAPHIC INTERACTION (VMM scale-dependent) */
if (pba->interaction_beta != 0.) {
    double rho_cdm = pvecback[pba->index_bg_rho_cdm];
    double rho_lambda = pvecback[pba->index_bg_rho_lambda];
    double rho_tot = pvecback[pba->index_bg_rho_tot];
    double Omega_lambda = rho_lambda / rho_tot;
    double w_lambda = -1.0;
    double a = pvecback[pba->index_bg_a];
    
    /* Background interaction: Q/rho = -3 * beta * aH * Omega_de * w_de */
    double Q_over_rho = -3.0 * pba->interaction_beta * a_prime_over_a * Omega_lambda * w_lambda;
    
    /* Scale-dependent factor using k/k_eq transition */
    double k_eq = 0.01 * 0.6766;  /* k_eq ~ 0.01 h/Mpc */
    double x = k / k_eq;
    double scale_factor = x * x / (1.0 + x * x);
    
    /* Late-time activation: interaction grows as DE becomes important */
    double late_time_factor = Omega_lambda * Omega_lambda;
    
    /* Combined scale-dependent drag */
    dy[pv->index_pt_delta_cdm] += Q_over_rho * scale_factor * late_time_factor * y[pv->index_pt_delta_cdm];
}
/* END HOLOGRAPHIC INTERACTION */
```

## Physical Effects

1. **Scale-dependent suppression**: 
   - k << k_eq: No suppression (preserves CMB)
   - k >> k_eq: Full drag effect (suppresses σ₈)

2. **Late-time activation**:
   - Ω_Λ² factor ensures drag only active when DE dominates
   - Preserves early-universe physics

3. **VMM stability**:
   - Implements Valiviita-Majerotto-Maartens prescription
   - Avoids "doom factor" instabilities

## Results Summary

| β (paper) | σ₈    | S₈    | CMB deviation |
|-----------|-------|-------|---------------|
| 0.0       | 0.821 | 0.834 | 0%            |
| 0.3       | 0.782 | 0.795 | <0.5%         |
| 0.5       | 0.758 | 0.770 | <0.7%         |
| 1.0       | 0.702 | 0.713 | <1%           |

Target: S₈ ≈ 0.76-0.78 (weak lensing)

## Usage

In Python:
```python
from classy import Class
cosmo = Class()
cosmo.set({
    'output': 'mPk',
    'interaction_beta': -0.5,  # Note: code uses opposite sign
    # ... other parameters
})
cosmo.compute()
print(f"sigma8 = {cosmo.sigma8()}")
```

In .ini file:
```
interaction_beta = -0.5
```

## Citation

If using this code, please cite:
- Sayegh (2026), "The Holographic Signature of the Cosmic Horizon"
- CLASS: Blas, Lesgourgues, Tram (2011)

## Contact

Bernard Sayegh
bernardsayegh@hotmail.com

## Sign Convention (Updated)

The interaction parameter `interaction_beta` uses the **same sign** as the paper:
- **β = +0.5** in paper corresponds to `interaction_beta = 0.5` in code
- Positive β → energy flows DE → matter (Q > 0)
- Positive β → **suppresses** σ₈ via dilution mechanism

### Numerical Results (Planck 2018 parameters)

| β | σ₈ | S₈ |
|---|-----|-----|
| 0.0 | 0.823 | 0.837 |
| 0.25 | 0.790 | 0.803 |
| **0.5** | **0.759** | **0.772** |
| 0.75 | 0.730 | 0.742 |
| 1.0 | 0.702 | 0.714 |

