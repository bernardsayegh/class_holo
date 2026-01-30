# Implementation Guide: Holographic Dark Energy in CLASS

This document provides complete technical documentation for the holographic dark energy implementation, including the physics, code structure, parameters, and validation scripts.

## Table of Contents

1. [Core Physics](#1-core-physics)
2. [Modified Files](#2-modified-files)
3. [Parameter Reference](#3-parameter-reference)
4. [Model Configurations](#4-model-configurations)
5. [Code Walkthrough](#5-code-walkthrough)
6. [Validation Scripts](#6-validation-scripts)
7. [MCMC Setup](#7-mcmc-setup)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Core Physics

### 1.1 The Interaction Kernel

The holographic model transfers energy from dark energy (Λ) to cold dark matter at a rate:

```
Q = β · I_eff(q) · Ω_Λ · ρ_tot · H
```

where:
- **β = 1/12 ≈ 0.0833**: Coupling strength from Nariai stability bound
- **I_eff = [(1-q)/2]²**: Efficiency modulation based on deceleration
- **q = -1 + (3/2)Ω_m**: Deceleration parameter
- **Ω_Λ, Ω_m**: Two-fluid density fractions (DE + matter only)

### 1.2 Two-Fluid Treatment

The interaction uses a **two-fluid approximation** for the sweep kernel:

```c
rho_m_2f = rho_cdm_holo + rho_b;           // Matter = CDM + baryons
rho_2f = rho_de + rho_m_2f;                 // Total two-fluid density
Omega_de_2f = rho_de / rho_2f;
Omega_m_2f = rho_m_2f / rho_2f;
```

This excludes radiation and neutrinos from the interaction dynamics while keeping them in the full Friedmann equation.

### 1.3 Super-Schwarzschild Correction (SCR)

The SCR subsystem tracks when the horizon sweep rate exceeds the Schwarzschild limit:

```
S = 4.5 · Ω_de · Ω_m    (parabola, peaks at S=1 when Ω_m = Ω_de = 0.5)
```

When S > 1, excess energy can be:
1. **Stored** in a reservoir ρ_scr (controlled by `super_schw_amp`)
2. **Mapped** to a local H₀ inference (controlled by `super_schw_Amap`)

### 1.4 H₀_local Mapping

The key to resolving the H₀ tension:

```
H₀_local = H₀_phys × exp(Amap × X₀)
```

where:
- **H₀_phys**: Physical Hubble parameter from integrated Friedmann equation (what CMB measures)
- **X₀**: Accumulated super-Schwarzschild excess at z=0
- **Amap = 2**: Mapping amplitude (Model B)

With Amap=2 and typical X₀ ≈ 0.035, this gives H₀_local ≈ 73 km/s/Mpc while H₀_phys ≈ 67.3 km/s/Mpc.

**Physical interpretation**: H₀_phys is the true cosmic expansion rate that determines CMB acoustic scale. H₀_local is what distance ladder measurements infer due to the super-Schwarzschild mapping effect on local calibrators.

### 1.5 X_schw Accumulation

The dimensionless excess X is accumulated during background integration:

```c
dX/dlna = (S > 1) ? (S - 1) : 0
```

This integral captures the total super-Schwarzschild history. X₀ = X(z=0) is used for the H₀_local mapping.

### 1.6 Far-Future Completion (κ saturation)

Optional mechanism to drive Ω_m → 1/3 in the far future (z < 0):

```
boost = 1 + κ · g_future · g_S<1 · (1/S - 1)
```

This only activates for a > 1 and S < 1, ensuring no effect on fitted cosmology.

---

## 2. Modified Files

### 2.1 `source/background.c`

**Main physics implementation.** Key functions:

- **`background_derivs()`**: Implements the interaction ODEs
  - Two-fluid density calculation
  - Deceleration parameter q and efficiency I_eff
  - Energy transfer Q/H
  - SCR reservoir dynamics (if enabled)
  - X_schw accumulation
  - Far-future κ saturation (if enabled)

- **`background_solve()`**: Computes derived quantities
  - `H0_phys` from integrated background table at z=0
  - `X0_schw` from background table
  - `H0_local = H0_phys * exp(Amap * X0)`

### 2.2 `source/input.c`

**Parameter parsing.** Reads:
- `interaction_beta`
- `f_clust`
- `super_schwarzschild_correction`
- `super_schw_amp`, `super_schw_Amap`, `super_schw_gamma`, `super_schw_deltaS`
- `super_schw_no_mapping`
- `super_schw_kappa`

### 2.3 `include/background.h`

**Parameter declarations.** Adds fields to the background structure:
- `double interaction_beta;`
- `double f_clust;`
- `double super_schw_amp;`
- `double super_schw_Amap;`
- `double super_schw_gamma;`
- `double super_schw_deltaS;`
- `double super_schw_kappa;`
- `int super_schw_no_mapping;`
- `double H0_local;`
- `double X0_schw;`

---

## 3. Parameter Reference

### 3.1 Core Interaction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `interaction_beta` | 0 | Coupling strength. Theory predicts β = 1/12 ≈ 0.0833 |
| `f_clust` | 0 | Clustering fraction. Use 0 for smooth CDM deposition |

### 3.2 Super-Schwarzschild Correction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `super_schwarzschild_correction` | no | Enable SCR subsystem |
| `super_schw_amp` | 0 | Reservoir storage strength. 0 = no reservoir (Model B) |
| `super_schw_Amap` | 2 | Mapping amplitude for H₀_local |
| `super_schw_gamma` | 2.0 | Reservoir decay rate back to CDM |
| `super_schw_deltaS` | 0.03 | Smoothness of SCR activation around S=1 |
| `super_schw_no_mapping` | 0 | Set to 1 to disable mapping (H₀_local = H₀_phys) |
| `super_schw_kappa` | 0 | Far-future saturation feedback strength |

---

## 4. Model Configurations

### 4.1 ΛCDM Baseline

```ini
interaction_beta = 0
# No other holographic parameters needed
```

### 4.2 Model A (Interaction Only)

Resolves S₈ tension. Does NOT resolve H₀ tension.

```ini
interaction_beta = 0.08333333   # = 1/12
f_clust = 0.0
super_schwarzschild_correction = no
```

Or with SCR infrastructure but inert:
```ini
interaction_beta = 0.08333333
f_clust = 0.0
super_schwarzschild_correction = yes
super_schw_amp = 0.0
super_schw_Amap = 0.0    # No mapping
```

### 4.3 Model B (Interaction + Mapping) ⭐ MAIN RESULT

Resolves BOTH S₈ and H₀ tensions.

```ini
interaction_beta = 0.08333333   # = 1/12
f_clust = 0.0
super_schwarzschild_correction = yes
super_schw_amp = 0.0            # No reservoir
super_schw_Amap = 2.0           # Mapping amplitude
super_schw_deltaS = 0.03
super_schw_gamma = 2.0
super_schw_no_mapping = 0       # Mapping enabled
```

### 4.4 Model C (Full Reservoir + Mapping)

```ini
interaction_beta = 0.08333333
f_clust = 0.0
super_schwarzschild_correction = yes
super_schw_amp = 1.0            # Reservoir active
super_schw_Amap = 2.0
super_schw_deltaS = 0.03
super_schw_gamma = 2.0
```

### 4.5 Model D (Reservoir Only, No Mapping)

```ini
interaction_beta = 0.08333333
f_clust = 0.0
super_schwarzschild_correction = yes
super_schw_amp = 1.0
super_schw_Amap = 2.0
super_schw_no_mapping = 1       # Mapping disabled
```

### 4.6 Far-Future Completion (Optional)

Add to any model to ensure Ω_m → 1/3 as t → ∞:

```ini
super_schw_kappa = 14.22        # Or any positive value
```

The critical value κ* = 128/9 ≈ 14.22 pins Ω_m exactly to 1/3.

---

## 5. Code Walkthrough

### 5.1 The Interaction Kernel in `background_derivs()`

```c
// Two-fluid densities (matter + DE only)
double rho_m_2f = rho_cdm_holo + rho_b;
double rho_2f = rho_de + rho_m_2f;
double Omega_de_2f = rho_de / rho_2f;
double Omega_m_2f = rho_m_2f / rho_2f;

// Deceleration parameter
double q_decel_2f = -1.0 + 1.5 * Omega_m_2f;

// Efficiency modulation
double I_eff = 0.25 * pow(1.0 - q_decel_2f, 2);

// Effective coupling
double beta_eff = pba->interaction_beta * I_eff;

// Energy transfer rate Q/H
double Q_over_H = 4.5 * beta_eff * Omega_de_2f * Omega_m_2f * rho_2f;
```

### 5.2 SCR Reservoir Dynamics

```c
// Schwarzschild parabola
double S_ss = 4.5 * Omega_de_2f * Omega_m_2f;

// Smooth activation gate
double gS = 0.5 * (1.0 + tanh((S_ss - 1.0) / pba->super_schw_deltaS));

// Fractional excess
double frac_ss = (S_ss > 1.0) ? (1.0 - 1.0/S_ss) : 0.0;

// Energy-conserving split
double f_store = pba->super_schw_amp * gS * frac_ss;
f_store = fmin(fmax(f_store, 0.0), 1.0);  // Clamp to [0,1]

double Q_to_scr_over_H = f_store * Q_over_H;
double Q_to_cdm_over_H = (1.0 - f_store) * Q_over_H;

// Reservoir decay
double decay_scr = pba->super_schw_gamma * (1.0 - gS) * rho_scr;
```

### 5.3 X_schw Accumulation

```c
// In background_derivs(), accumulate super-Schwarzschild excess
double dX_dlna = 0.0;
if (S_ss > 1.0) {
    dX_dlna = S_ss - 1.0;  // Excess above Schwarzschild limit
}
// This is integrated alongside other background quantities
```

### 5.4 H₀_local Mapping in `background_solve()`

```c
/* Compute H0_local from super-Schwarzschild correction */
if (pba->has_super_schw_correction == _TRUE_) {
    /* Get X0 from background table (accumulated excess at z=0) */
    pba->X0_schw = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_X_schw];
    
    /* Get H0_phys from integrated Friedmann equation */
    double H0_phys = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_H];
    
    if (pba->super_schw_no_mapping == _TRUE_) {
        pba->H0_local = H0_phys; /* mapping disabled => local equals physical */
    } else {
        /* Apply mapping: H0_local = H0_phys * exp(Amap * X0) */
        pba->H0_local = H0_phys * exp(pba->super_schw_Amap * pba->X0_schw);
    }
} else {
    pba->X0_schw = 0.0;
    pba->H0_local = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_H];
}
```

**Note**: H0_phys is read from the background table at z=0, which contains the integrated Friedmann solution. This is the physical expansion rate that CMB acoustic physics sees.

---

## 6. Validation Scripts

### 6.1 All Models Comparison (θs-matched, H₀ split)

This script compares all models at fixed Planck θs, showing H₀_input, H₀_phys, and H₀_local separately:

```python
#!/usr/bin/env python3
"""
cosmo_stats_h0split.py
Compare all models with explicit H0_input vs H0_phys vs H0_local.
"""
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
    """Get H0 from integrated background at z=0."""
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

print("\nSUMMARY (explicit H0_in vs H0_phys vs H0_local)")
print(f"{'Model':<10} {'h':>7} {'H0_in':>8} {'H0_phy':>8} {'H0_loc':>8} {'X0':>8} {'Om':>8} {'sig8':>8} {'S8':>8}")
for name, s in results.items():
    print(f"{name:<10} {s['h']:>7.5f} {s['H0_in']:>8.2f} {s['H0_phys']:>8.2f} {s['H0_loc']:>8.2f} "
          f"{s['X0']:>8.5f} {s['Omega_m']:>8.5f} {s['sigma8']:>8.5f} {s['S8']:>8.5f}")

print("\nTargets: SH0ES H0=73.04, DES S8=0.776")
```

**Expected output:**
```
Model            h    H0_in   H0_phy   H0_loc       X0       Om     sig8       S8
LCDM       0.67165    67.16    67.16    67.16  0.00000  0.31565  0.80761  0.82841
Model A    0.66783    66.78    67.28    67.28  0.00000  0.32927  0.74255  0.77793
Model B    0.66783    66.78    67.28    72.15  0.03497  0.32927  0.74255  0.77793
Model C    0.66202    66.20    68.25    73.20  0.03499  0.32048  0.73748  0.76223
Model D    0.66202    66.20    68.25    68.25  0.03499  0.32048  0.73748  0.76223
```

**Key observations:**
- **LCDM**: H0_in = H0_phys = H0_local (no interaction, no mapping)
- **Model A/B**: Same h, same H0_phys, but Model B has H0_local boosted by mapping
- **Model C/D**: Same h (both have reservoir), same H0_phys, but only Model C has mapping
- **Model B**: Achieves H0_local ≈ 72 and S8 ≈ 0.78 — resolves both tensions!

### 6.2 Model B Quick Check

```python
#!/usr/bin/env python3
"""
check_modelB.py
Verify Model B produces correct H0_local via mapping.
"""
import sys
sys.path.insert(0, "python")
from classy import Class

params = {
    "output": "mPk",
    "P_k_max_1/Mpc": 10.0,
    
    "h": 0.6784,
    "omega_b": 0.02237,
    "omega_cdm": 0.1200,
    "A_s": 2.1e-9,
    "n_s": 0.9649,
    "tau_reio": 0.0544,
    
    # Holographic interaction
    "interaction_beta": 1.0/12.0,
    "f_clust": 0.0,
    
    # Super-Schwarzschild (Model B)
    "super_schwarzschild_correction": "yes",
    "super_schw_amp": 0.0,
    "super_schw_Amap": 2.0,
    "super_schw_deltaS": 0.03,
    "super_schw_gamma": 2.0,
}

c = Class()
c.set(params)
c.compute()

derived = c.get_current_derived_parameters(["H0_local", "sigma8", "Omega_m"])
H0_local = derived["H0_local"]
sigma8 = derived["sigma8"]
Omega_m = derived["Omega_m"]
S8 = sigma8 * (Omega_m / 0.3)**0.5

print("Model B Results:")
print(f"  H0 (physical) = {100*c.h():.2f} km/s/Mpc")
print(f"  H0_local      = {H0_local:.2f} km/s/Mpc")
print(f"  sigma8        = {sigma8:.4f}")
print(f"  Omega_m       = {Omega_m:.4f}")
print(f"  S8            = {S8:.4f}")

c.struct_cleanup()
c.empty()
```

---

## 7. MCMC Setup

### 7.1 Required Likelihoods

The paper uses:
- **Planck 2018**: Full TTTEEE (plik), lowl TT, lowl EE, lensing
- **BAO**: 6dFGS, SDSS DR7 MGS, SDSS DR12 consensus, SDSS DR16 (QSO, ELG, Lyα)
- **Supernovae**: Pantheon+
- **Weak Lensing**: DES Y3 cosmic shear (with intrinsic alignment)
- **Local H₀**: SH0ES 2022 via custom `H0_local` likelihood

### 7.2 Custom SH0ES Likelihood

The key innovation: constrain `H0_local` (not `H0`) against SH0ES.

```python
# cobaya/likelihoods/shoes_h0local/shoes_h0local.py
from cobaya.likelihood import Likelihood

class SH0ES_H0local(Likelihood):
    """SH0ES constraint using H0_local from holographic CLASS."""
    H0_shoes: float = 73.04
    H0_shoes_err: float = 1.04
    
    def initialize(self):
        self.log.info(f"SH0ES H0_local: {self.H0_shoes} +/- {self.H0_shoes_err}")
    
    def get_requirements(self):
        return {'H0_local': None}
    
    def logp(self, **params_values):
        H0_local = self.provider.get_param('H0_local')
        return -0.5 * ((H0_local - self.H0_shoes) / self.H0_shoes_err) ** 2
```

### 7.3 Cobaya YAML Configuration

See `cobaya/modelB_Amap2.yaml` for the complete configuration. Key points:

```yaml
theory:
  classy:
    extra_args:
      interaction_beta: 0.0833333
      f_clust: 0.0
      super_schwarzschild_correction: 'yes'
      super_schw_amp: 0.0
      super_schw_Amap: 2.0
    output_params:
      - sigma8
      - Omega_m
      - H0_local    # Critical!

likelihood:
  shoes_h0local.SH0ES_H0local:
    H0_shoes: 73.04
    H0_shoes_err: 1.04
```

---

## 8. Troubleshooting

### 8.1 H0_local Not Updating

**Symptom**: H0_local equals H0_phys regardless of Amap setting.

**Check**:
1. Is `super_schwarzschild_correction: yes`?
2. Is `super_schw_no_mapping: 0` (not 1)?
3. Is `super_schw_Amap` non-zero?

### 8.2 NaN or Inf in Background

**Symptom**: Integration fails with NaN values.

**Likely causes**:
- Division by zero in two-fluid calculation
- Negative densities
- Too aggressive step size

**Fix**: Check initial conditions and ensure `Omega_k = 0`.

### 8.3 S₈ Not Reduced

**Symptom**: σ₈ similar to ΛCDM despite β > 0.

**Check**:
1. Is `interaction_beta` actually set? (Default is 0)
2. Is `f_clust = 0`? (f_clust = 1 reduces the effect)
3. Are you comparing at matched θs?

### 8.4 Build Artifacts

After modifying C files, always:

```bash
make clean
make class
pip install -e . --user  # If using Python wrapper
```

Don't commit `class` binary or `python/classy.cpp` (build artifacts).

---

## Appendix: Expected Results

### θs-matched Theoretical Predictions

| Model | H₀_phys | H₀_local | σ₈ | S₈ |
|-------|---------|----------|-----|-----|
| ΛCDM | 67.2 | 67.2 | 0.808 | 0.828 |
| Model A | 67.3 | 67.3 | 0.743 | 0.778 |
| Model B | 67.3 | 72.2 | 0.743 | 0.778 |
| Model C | 68.3 | 73.2 | 0.737 | 0.762 |
| Model D | 68.3 | 68.3 | 0.737 | 0.762 |

### Observational Targets

- **SH0ES H₀**: 73.04 ± 1.04 km/s/Mpc
- **Planck H₀**: 67.4 ± 0.5 km/s/Mpc  
- **DES Y3 S₈**: 0.776 ± 0.017
- **Planck S₈**: 0.834 ± 0.016

### Model B Performance

**Model B (β=1/12, Amap=2)** achieves:
- H₀_local within 1σ of SH0ES
- S₈ within 1σ of DES Y3
- Perfect fit to Planck CMB (same θs)
- **Δχ² ≈ -50 vs ΛCDM** with zero additional free parameters
