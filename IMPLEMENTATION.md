# Implementation Guide: Holographic Dark Energy in CLASS

Complete technical documentation for the holographic dark energy implementation in CLASS, including physics, code structure, parameters, validation, and MCMC setup.

## Table of Contents

1. [Core Physics](#1-core-physics)
2. [Modified Files](#2-modified-files)
3. [Parameter Reference](#3-parameter-reference)
4. [Model Configurations](#4-model-configurations)
5. [Code Walkthrough](#5-code-walkthrough)
6. [Validation Scripts](#6-validation-scripts)
7. [MCMC Setup](#7-mcmc-setup)
8. [Troubleshooting](#8-troubleshooting)
9. [Appendix: Diagnostic Results](#9-appendix-diagnostic-results)

---

## 1. Core Physics

### 1.1 The Interaction Kernel

The holographic model transfers energy from dark energy (Λ) to cold dark matter at a rate:

```
Q = β · I_eff(q) · Ω_Λ · ρ_tot · H
```

where:
- **β = 1/12 ≈ 0.0833**: Coupling strength from the Nariai stability bound
- **I_eff**: Efficiency modulation based on deceleration parameter (see §1.2)
- **q**: Deceleration parameter
- **Ω_Λ, Ω_m**: Two-fluid density fractions (DE + matter only, excluding radiation/neutrinos)

### 1.2 I_eff Modulation Types

The efficiency factor I_eff controls how strongly the interaction couples at different epochs. Five formulations are implemented, selected by `interaction_ieff_type`:

| ieff_type | Name | Formula | Attractor value (Ω_Λ=2/3) | Notes |
|-----------|------|---------|---------------------------|-------|
| 0 | Thermodynamic squared | [(1−q)/2]² | 0.5625 | Original formulation |
| 1 | Lagrangian squared | [3/(2−q)]² | 0.36 | Over-suppresses; needs β ≈ 1/72 |
| 2 | Lagrangian linear | 3/(2−q) | 0.6 | Over-suppresses; needs β ≈ 1/36 |
| 3 | Thermodynamic fourth power | [(1−q)/2]⁴ | 0.3164 | Too weak; barely affects cosmology |
| **4** | **Thermodynamic linear (unsquared)** | **(1−q)/2** | **0.75** | **Preferred: simplest, best fits** |

The two-fluid deceleration parameter is:
```
q_2f = -1 + (3/2) · Ω_m,2f
```

where Ω_m,2f = ρ_m / (ρ_m + ρ_de), excluding radiation and neutrinos.

**Diagnostic results at fixed β=1/12, Amap=2:**

| ieff_type | σ₈ | S₈ | Ω_m | H₀ |
|-----------|------|------|------|------|
| 0 (squared) | 0.7884 | 0.8192 | 0.3239 | 67.87 |
| 4 (unsquared) | 0.7666 | 0.8051 | 0.3309 | 68.22 |
| 1 (lag²) | 0.6309 | 0.7116 | 0.3817 | 70.97 |
| 2 (lag) | 0.6873 | 0.7518 | 0.3589 | 69.69 |

The unsquared (ieff_type=4) is the preferred formulation: it is the simplest (linear surface gravity ratio), produces stronger S₈ suppression than squared without over-suppressing like Lagrangian variants, and achieves the best MCMC fits (Δχ² ≈ −48 vs ΛCDM).

### 1.3 Two-Fluid Treatment

The interaction uses a **two-fluid approximation** for the sweep kernel:

```c
rho_m_2f = rho_cdm_holo + rho_b;        // Matter = CDM + baryons
rho_2f   = rho_de + rho_m_2f;            // Total two-fluid density
Omega_de_2f = rho_de / rho_2f;
Omega_m_2f  = rho_m_2f / rho_2f;
```

This excludes radiation and neutrinos from the interaction dynamics while keeping them in the full Friedmann equation. Critical: earlier versions incorrectly used total-density fractions including radiation, which caused significant overestimation of perturbation effects. The two-fluid treatment must be used consistently in both background.c and perturbations.c.

### 1.4 Super-Schwarzschild Correction (SCR)

The SCR subsystem tracks when the horizon sweep rate exceeds the Schwarzschild limit:

```
S = 4.5 · Ω_de · Ω_m    (parabola, peaks at S=1 when Ω_m = Ω_de = 0.5)
```

When S > 1, excess energy can be:
1. **Stored** in a reservoir ρ_scr (controlled by `super_schw_amp`)
2. **Mapped** to a local H₀ inference (controlled by `super_schw_Amap`)

### 1.5 X_schw: ODE Formulations

The dimensionless super-Schwarzschild excess X determines the H₀ mapping strength. Three formulations are implemented, selected by `super_schw_ode`:

**ODE with dilution (super_schw_ode = 1):**
```
dΨ/d(ln a) = −(1+q)(Ψ + Δ)
where Δ = max(S − 1, 0)
```
Standard perturbation equation derived from perturbing the Misner-Sharp energy: define δE = −ΨE⁰, take the product rule, use d(ln E⁰)/d(ln a) = 1+q. The dilution term (1+q)Ψ arises because the background horizon energy grows — a fixed absolute perturbation becomes a smaller fractional departure. Gives |Ψ| ≈ 0.022 at z=0.

**ODE without dilution (super_schw_ode = 2):**
```
dΨ/d(ln a) = −(1+q)Δ
```
Retains the (1+q) source weighting (horizon growth rate determines how much unprocessed flux is produced per e-fold) but drops the dilution term. Physically: the source *should* scale with horizon growth, but the accumulated record doesn't decay because once entropy is produced it is permanent (second law). The denominator against which Ψ is normalised is fixed at the time of production, not updating with E⁰. Gives |Ψ| ≈ 0.027 at z=0 — intermediate between full ODE and accumulator.

**Accumulator (super_schw_ode = 0, legacy):**
```
dX/d(ln a) = Δ
```
Pure monotonic integral of the excess — no (1+q) weighting, no dilution. X₀ ≈ 0.035 at z=0. Treats X as a raw entropy ledger rather than a geometric perturbation. Retained for backward compatibility and comparison with earlier results.

**Diagnostic comparison (β=1/12, unsquared, Amap=2):**

| Mode | X₀ | H₀_loc (Amap=2) |
|------|------|-----------------|
| Accumulator (ode=0) | 0.0355 | ~72.6 |
| ODE no-decay (ode=2) | 0.0266 | ~71.6 |
| ODE with dilution (ode=1) | 0.0225 | ~71.1 |

σ₈ and S₈ are identical across all three modes — X₀ only affects the mapping, not perturbations.

**Theoretical comparison:**
- **ODE with dilution**: Most conservative. Every step is standard perturbation theory. Ψ is a fractional perturbation δE/E⁰ and the denominator updates as E⁰ grows. Easiest to justify to referees.
- **ODE without dilution**: Intermediate. Source weighted by horizon growth rate (physical), but accumulated record is permanent (entropy monotonicity). Ψ tracks cumulative entropy production whose geometric significance is fixed at production time.
- **Accumulator**: Strongest mapping. Treats X as a pure bookkeeping quantity with no geometric weighting. Requires additional physical arguments (path-integral distance ladder, separate-universe comparison).

### 1.6 H₀_local Mapping

The key to resolving the H₀ tension:

```
H₀_local = H₀_phys × exp(Amap × X₀)
```

where:
- **H₀_phys**: Physical Hubble parameter from the integrated Friedmann equation (what CMB measures)
- **X₀**: Super-Schwarzschild excess at z=0 (from ODE with dilution, ODE no-decay, or accumulator)
- **Amap**: Mapping amplitude (theory predicts Amap=2)

With Amap=2, unsquared I_eff, and β=1/12:
- ODE with dilution: X₀ ≈ 0.022 → H₀_local ≈ 71.1 km/s/Mpc
- ODE no-decay: X₀ ≈ 0.027 → H₀_local ≈ 71.6 km/s/Mpc
- Accumulator: X₀ ≈ 0.035 → H₀_local ≈ 72.6 km/s/Mpc

**Physical interpretation**: H₀_phys is the true cosmic expansion rate that determines the CMB acoustic scale. H₀_local is what distance ladder measurements infer due to the super-Schwarzschild mapping effect on local calibrators.

### 1.7 Clustering Fraction (f_clust)

Controls how the injected dark energy perturbation is distributed:

```
δρ_de_injected = f_clust × (clustered) + (1−f_clust) × (homogeneous)
```

- **f_clust = 0** (default, theory prediction): Injection is homogeneous — the horizon cannot resolve sub-horizon structure
- **f_clust = 1**: Injection clusters with matter — kills S₈ suppression

The holographic principle predicts f_clust = 0 from first principles. MCMC validation:
- f_clust = 1 fixed: Δχ² = +146 vs ΛCDM (catastrophic)
- f_clust free: converges to f_clust ≈ 0.03 ± 0.01 (consistent with zero)

### 1.8 CPL Surrogate Analysis

The holographic background is indistinguishable from ΛCDM. Fitting H(z) and D_M(z) to CPL parameterisation w(a) = w₀ + wₐ(1−a):

| I_eff type | w₀ | wₐ | Ω_m,CPL | Hidden matter | Fit quality |
|------------|------|------|---------|---------------|-------------|
| Squared | −1.004 | −0.008 | 0.315 | 4.4% | 0.008% |
| Unsquared | −0.995 | −0.023 | 0.316 | 7.2% | 0.004% |

The CPL surrogate produces σ₈ and S₈ identical to ΛCDM (δS₈ = 0.0000). The entire S₈ suppression (δS₈ ≈ −0.02 to −0.03) comes purely from the perturbation-level interaction. An observer fitting only background expansion would infer standard ΛCDM, yet σ₈ is 3–4% lower.

---

## 2. Modified Files

### 2.1 `source/background.c`

Main modifications:
- **Two-fluid density computation**: Computes Ω_de,2f and Ω_m,2f excluding radiation/neutrinos
- **I_eff switch**: Implements all five ieff_type variants (0–4)
- **Interaction rate Q**: Computes Q = β · I_eff · Ω_Λ · ρ_tot · H
- **SCR parabola**: S = 4.5 · Ω_de · Ω_m
- **X_schw integration**: ODE with dilution (ode=1), ODE no-decay (ode=2), and accumulator (ode=0, legacy)
- **Reservoir ρ_scr**: Optional energy storage when S > 1
- **H₀_local derived parameter**: H₀_phys × exp(Amap × X₀)
- **Background table columns**: X_schw, Q_interaction, I_eff, S_parabola

### 2.2 `source/perturbations.c`

Main modifications:
- **Perturbation coupling Q_pert**: Uses same two-fluid fractions and I_eff as background
- **δρ_de source term**: Homogeneous injection (f_clust=0) or clustered (f_clust=1)
- **Consistent I_eff**: Same ieff_type switch as background.c — must match exactly
- **Order of execution**: Q_scr_to_cdm computed before use in density/velocity equations

### 2.3 `include/background.h`

Added fields:
- `interaction_beta`, `f_clust`, `interaction_ieff_type`
- `super_schwarzschild_correction`, `super_schw_amp`, `super_schw_Amap`
- `super_schw_gamma`, `super_schw_deltaS`, `super_schw_kappa`
- `super_schw_ode`
- Background table indices for new columns

### 2.4 `source/input.c`

Reads all new parameters from .ini/.yaml files with defaults.

---

## 3. Parameter Reference

### 3.1 Interaction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `interaction_beta` | 0 | Coupling strength. Theory predicts 1/12 ≈ 0.0833. Set to 0 to disable interaction. |
| `interaction_ieff_type` | 0 | I_eff formulation: 0=squared, 1=lag², 2=lag, 3=fourth power, **4=unsquared (preferred)** |
| `f_clust` | 0 | Clustering fraction [0,1]. Theory predicts 0 (homogeneous injection). |

### 3.2 Super-Schwarzschild Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `super_schwarzschild_correction` | no | Master switch: 'yes' to enable SCR subsystem |
| `super_schw_amp` | 0 | Reservoir amplitude. 0=no reservoir (Models A,B). 1=reservoir active (Models C,D). |
| `super_schw_Amap` | 2 | H₀ mapping amplitude. 0=no mapping (Models A,D). 2=theory prediction (Models B,C). |
| `super_schw_gamma` | 2.0 | SCR parabola exponent |
| `super_schw_deltaS` | 0.03 | SCR transition width |
| `super_schw_ode` | 0 | X formulation: 0=accumulator (legacy), **1=ODE with dilution**, **2=ODE without dilution** |
| `super_schw_kappa` | 0 | Far-future saturation (set to 0) |

### 3.3 Derived Output Parameters

| Parameter | Description |
|-----------|-------------|
| `H0_local` | H₀_phys × exp(Amap × X₀). This is what SH0ES measures. |
| `X0_schw` | Accumulated super-Schwarzschild excess at z=0 |
| `sigma8` | RMS matter fluctuations in 8 Mpc/h spheres |
| `Omega_m` | Total matter density fraction |
| `S8` | σ₈ × √(Ω_m/0.3) — the weak lensing observable |

---

## 4. Model Configurations

### 4.1 ΛCDM Baseline

```yaml
extra_args:
  interaction_beta: 0
```

No interaction. Standard cosmology. All SCR parameters irrelevant.

### 4.2 Model A (Perturbation Damping Only)

Resolves S₈ only. No H₀ mapping, no reservoir.

```yaml
extra_args:
  interaction_beta: 0.0833
  interaction_ieff_type: 4
  f_clust: 0.0
```

Does NOT resolve H₀. H₀_local = H₀_phys (no mapping).

### 4.3 Model B (Damping + H₀ Mapping) ⭐ MAIN RESULT

Resolves BOTH S₈ and H₀ tensions.

```yaml
extra_args:
  interaction_beta: 0.0833
  interaction_ieff_type: 4
  f_clust: 0.0
  super_schwarzschild_correction: 'yes'
  super_schw_amp: 0.0
  super_schw_Amap: 2.0
  super_schw_gamma: 2.0
  super_schw_deltaS: 0.03
```

For ODE with dilution, add:
```yaml
  super_schw_ode: 1
```

For ODE without dilution (no-decay), add:
```yaml
  super_schw_ode: 2
```

### 4.4 Model C (Damping + Reservoir + Mapping)

```yaml
extra_args:
  interaction_beta: 0.0833
  interaction_ieff_type: 4
  f_clust: 0.0
  super_schwarzschild_correction: 'yes'
  super_schw_amp: 1.0          # Reservoir active
  super_schw_Amap: 2.0          # Mapping active
  super_schw_gamma: 2.0
  super_schw_deltaS: 0.03
```

### 4.5 Model D (Reservoir Only, No Mapping)

```yaml
extra_args:
  interaction_beta: 0.0833
  interaction_ieff_type: 4
  f_clust: 0.0
  super_schwarzschild_correction: 'yes'
  super_schw_amp: 1.0          # Reservoir active
  super_schw_Amap: 0.0          # No mapping
  super_schw_gamma: 2.0
  super_schw_deltaS: 0.03
```

### 4.6 acdm (Mapping-Only Control)

Maps H₀ without the dark energy interaction. Tests whether the mapping mechanism alone improves fits. Uses standard ΛCDM perturbations.

```yaml
extra_args:
  super_schwarzschild_correction: 'yes'
  super_schw_amp: 0.0
  super_schw_Amap: 2.0
  super_schw_gamma: 2.0
  super_schw_deltaS: 0.03
  super_schw_ode: 2             # Match ODE mode to Model B comparison
  super_schw_no_mapping: 0
```

Note: acdm without interaction_beta produces S₈ identical to ΛCDM (mapping is background-only, does not affect perturbations), but still gives H₀_local boost. The `super_schw_no_mapping: 0` flag ensures DES Y3 is evaluated in the local frame for frame-control runs.

### 4.7 Model Summary Table

| Model | interaction_beta | ieff_type | f_clust | SCR | amp | Amap | Resolves |
|-------|-----------------|-----------|---------|-----|-----|------|----------|
| ΛCDM | 0 | — | — | no | — | — | — |
| acdm | 0 | — | — | yes | 0 | 2 | H₀ only |
| A | 1/12 | 4 | 0 | no | — | 0 | S₈ only |
| **B** | **1/12** | **4** | **0** | **yes** | **0** | **2** | **S₈ + H₀** |
| C | 1/12 | 4 | 0 | yes | 1 | 2 | S₈ + H₀ + reservoir |
| D | 1/12 | 4 | 0 | yes | 1 | 0 | S₈ + reservoir |

---

## 5. Code Walkthrough

### 5.1 background.c — Key Sections

**Two-fluid densities (in background_functions):**
```c
rho_m_2f = pvecback[pba->index_bg_rho_cdm] + pvecback[pba->index_bg_rho_b];
rho_de   = pvecback[pba->index_bg_rho_lambda];
rho_2f   = rho_m_2f + rho_de;
Omega_m_2f  = rho_m_2f / rho_2f;
Omega_de_2f = rho_de / rho_2f;
```

**I_eff switch:**
```c
switch(pba->interaction_ieff_type) {
  case 0:  // Thermodynamic squared
    q_2f = -1.0 + 1.5 * Omega_m_2f;
    I_eff = 0.25 * (1.0 - q_2f) * (1.0 - q_2f);
    break;
  case 4:  // Thermodynamic linear (unsquared) — PREFERRED
    q_2f = -1.0 + 1.5 * Omega_m_2f;
    I_eff = 0.5 * (1.0 - q_2f);
    break;
  case 1:  // Lagrangian squared
    q_2f = -1.0 + 1.5 * Omega_m_2f;
    I_eff = 9.0 / ((2.0 - q_2f) * (2.0 - q_2f));
    break;
  case 2:  // Lagrangian linear
    q_2f = -1.0 + 1.5 * Omega_m_2f;
    I_eff = 3.0 / (2.0 - q_2f);
    break;
  case 3:  // Thermodynamic fourth power
    q_2f = -1.0 + 1.5 * Omega_m_2f;
    double half_1mq = 0.5 * (1.0 - q_2f);
    I_eff = half_1mq * half_1mq * half_1mq * half_1mq;
    break;
}
```

**X_schw ODE modes:**
```c
if (pba->super_schw_ode == 1) {
  // ODE with dilution: dPsi/dlna = -(1+q)*(Psi + Delta)
  double Delta = (S > 1.0) ? (S - 1.0) : 0.0;
  double qp1 = 1.0 + q_decel;
  dX_dlna = qp1 * Delta - qp1 * X_current;
} else if (pba->super_schw_ode == 2) {
  // ODE without dilution: dPsi/dlna = -(1+q)*Delta
  double Delta = (S > 1.0) ? (S - 1.0) : 0.0;
  double qp1 = 1.0 + q_decel;
  dX_dlna = qp1 * Delta;
} else {
  // Accumulator (legacy): dX/dlna = max(S-1, 0)
  dX_dlna = (S > 1.0) ? (S - 1.0) : 0.0;
}
```

### 5.2 perturbations.c — Key Sections

The perturbation coupling must use identical two-fluid fractions and I_eff as background:

```c
// Two-fluid fractions (NOT total density)
rho_m_2f = pvecback[pba->index_bg_rho_cdm] + pvecback[pba->index_bg_rho_b];
rho_de   = pvecback[pba->index_bg_rho_lambda];
rho_2f   = rho_m_2f + rho_de;

// Same I_eff switch as background.c
switch(pba->interaction_ieff_type) {
  case 0: I_eff = 0.25 * (1-q_2f)*(1-q_2f); break;
  case 4: I_eff = 0.5 * (1-q_2f); break;
  // ... etc
}

Q_pert = pba->interaction_beta * I_eff * Omega_de_2f * rho_2f * a_prime_over_a;
```

**Critical bug history**: An earlier version computed Q_pert using total-density fractions (including radiation), causing significant overestimation. The fix uses two-fluid fractions consistently. This bug was the reason for withdrawing the initial PRD submission.

### 5.3 Derived Parameters

```c
// In background_output_derived:
H0_local = H0_phys * exp(Amap * X0);
// Stored as derived parameter accessible to likelihoods
```

---

## 6. Validation Scripts

### 6.1 I_eff Diagnostic

Tests all I_eff formulations at fixed cosmology:

```bash
cd ~/class_holo_test && python3 << 'PY'
import sys; sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np

THETA_S_TARGET = 1.040423
c_km_s = 299792.458
BASE = {
    "output": "mPk", "P_k_max_1/Mpc": 10.0,
    "omega_b": 0.02242, "omega_cdm": 0.11933,
    "n_s": 0.9665, "A_s": 2.1e-9, "tau_reio": 0.054,
    "N_ur": 2.0328, "N_ncdm": 1, "m_ncdm": 0.06, "Omega_k": 0.0,
}

def run(label, beta, ieff):
    params = dict(BASE)
    params["interaction_beta"] = float(beta)
    params["interaction_ieff_type"] = int(ieff)
    params["f_clust"] = 0.0
    def res(h):
        p = dict(params); p["h"] = float(h)
        c = Class(); c.set(p); c.compute()
        th = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
        c.struct_cleanup(); c.empty()
        return th - THETA_S_TARGET
    h = brentq(res, 0.55, 0.80, xtol=1e-8, rtol=1e-12)
    params["h"] = float(h)
    c = Class(); c.set(params); c.compute()
    H0 = c.Hubble(0) * c_km_s
    Om = c.Omega_m()
    s8 = c.sigma8()
    S8 = s8 * np.sqrt(Om/0.3)
    print(f"  {label:<30s}  h={h:.5f}  H0={H0:.3f}  Om={Om:.4f}  s8={s8:.4f}  S8={S8:.4f}")
    c.struct_cleanup(); c.empty()

print("I_eff diagnostic (beta=1/12, no mapping)")
for ieff, name in [(0,"squared"), (4,"unsquared"), (1,"lag^2"), (2,"lag"), (3,"fourth")]:
    run(f"ieff={ieff} ({name})", 1/12, ieff)
PY
```

### 6.2 CPL Surrogate

Verifies that S₈ suppression is entirely from the interaction, not background modification. See `cpl_base.py` and `cpl_class_verify.py`.

### 6.3 Model B Quick Check (all ODE modes)

```bash
cd ~/class_holo_test && python3 << 'PY'
import sys; sys.path.insert(0, "python")
from classy import Class
import numpy as np
c_km_s = 299792.458
base = {
    "output": "mPk", "P_k_max_1/Mpc": 10.0,
    "h": 0.6678, "omega_b": 0.02237, "omega_cdm": 0.1200,
    "A_s": 2.1e-9, "n_s": 0.9649, "tau_reio": 0.0544,
    "interaction_beta": 1.0/12.0,
    "interaction_ieff_type": 4,
    "f_clust": 0.0,
    "super_schwarzschild_correction": "yes",
    "super_schw_amp": 0.0, "super_schw_Amap": 2.0,
    "super_schw_deltaS": 0.03, "super_schw_gamma": 2.0,
}
print(f"{'Mode':20s} {'X0':>8} {'H0_phys':>8} {'H0_loc':>8} {'sigma8':>8} {'S8':>8}")
for ode, label in [(0,"accumulator"), (1,"ODE dilution"), (2,"ODE no-decay")]:
    p = dict(base); p["super_schw_ode"] = ode
    c = Class(); c.set(p); c.compute()
    bg = c.get_background()
    X0 = bg['X_schw'][-1]
    H0 = 100*c.h()
    H0_loc = H0 * np.exp(2.0 * X0)
    s8 = c.sigma8(); Om = c.Omega_m()
    S8 = s8 * np.sqrt(Om/0.3)
    print(f"{label:20s} {X0:8.5f} {H0:8.3f} {H0_loc:8.3f} {s8:8.4f} {S8:8.4f}")
    c.struct_cleanup(); c.empty()
PY
```

---

## 7. MCMC Setup

### 7.1 Required Likelihoods

- **CMB**: Planck 2018 full (TTTEEE high-ℓ, low-ℓ TT, low-ℓ EE, lensing)
- **BAO**: 6dFGS, SDSS DR7 MGS, SDSS DR12 consensus, SDSS DR16 (QSO, ELG, Lyα)
- **Supernovae**: Pantheon+
- **Weak Lensing**: DES Y3 cosmic shear (with intrinsic alignment nuisance: A_IA, alpha_IA)
- **Local inference**: SH0ES 2022 via custom `H0_local` likelihood

### 7.2 Custom SH0ES Likelihood

Constrains `H0_local` (not `H0`) against SH0ES:

```python
# cobaya/likelihoods/shoes_h0local/shoes_h0local.py
from cobaya.likelihood import Likelihood

class SH0ES_H0local(Likelihood):
    """SH0ES constraint on H0_local from holographic model."""
    h0_shoes: float = 73.04
    h0_shoes_err: float = 1.04

    def initialize(self):
        self.log.info(f"SH0ES H0_local: {self.h0_shoes} +/- {self.h0_shoes_err}")

    def get_requirements(self):
        return {'H0_local': None}

    def logp(self, **params_values):
        H0_local = self.provider.get_param('H0_local')
        return -0.5 * ((H0_local - self.h0_shoes) / self.h0_shoes_err) ** 2
```

### 7.3 Cobaya YAML Template (Model B, unsquared, ODE)

```yaml
theory:
  classy:
    extra_args:
      interaction_beta: 0.0833
      interaction_ieff_type: 4
      f_clust: 0.0
      super_schwarzschild_correction: 'yes'
      super_schw_amp: 0.0
      super_schw_Amap: 2.0
      super_schw_gamma: 2.0
      super_schw_deltaS: 0.03
      super_schw_ode: 1          # 1=ODE with dilution, 2=ODE no-decay
    output_parameters:
      - sigma8
      - Omega_m
      - H0_local

likelihood:
  shoes_h0local.SH0ES_H0local:
    h0_shoes: 73.04
    h0_shoes_err: 1.04
  # ... plus Planck, BAO, Pantheon+, DES Y3
```

### 7.4 Sampled Parameters

Standard 6-parameter ΛCDM basis plus nuisance:
- `omega_b`, `omega_cdm`, `h`, `A_s` (or `logA`), `n_s`, `tau_reio`
- `A_IA`, `alpha_IA` (DES Y3 intrinsic alignment)

Optional additional sampled parameters:
- `interaction_beta`: Sample [0, 0.3] for β-free runs
- `super_schw_Amap`: Sample [0, 10] for Amap-free runs
- `super_schw_amp`: Sample [0, 5] for reservoir amplitude-free runs
- `f_clust`: Sample [0, 1] for clustering fraction test

### 7.5 Convergence

Target R−1 < 0.01 for publication. Use `cobaya-run -r` to resume chains. Monitor with `checkmc.py`.

---

## 8. Troubleshooting

### 8.1 H0_local Not Updating

**Symptom**: H0_local equals H0_phys (no mapping effect).

**Check**:
- Is `super_schwarzschild_correction: 'yes'`?
- Is `super_schw_kappa: 0` (not 1)?
- Is `super_schw_Amap` non-zero?

### 8.2 NaN or Inf in Background

**Symptom**: Integration failure with NaN values.

**Likely causes**:
- Division by zero in density fractions
- Negative density fractions
- Too aggressive step size

**Fix**: Check initial conditions and ensure `Omega_k = 0`.

### 8.3 S₈ Not Reduced

**Symptom**: σ₈ similar to ΛCDM despite β > 0.

**Check**:
- Is `interaction_beta` properly set? (default is 0)
- Is `f_clust = 0`? (f_clust = 1 kills suppression)
- Have you matched θ_s? (different h changes Ω_m which affects S₈)

### 8.4 Parameter Appears in Both Input and Extra Args

**Symptom**: Cobaya error about repeated parameter definitions.

**Fix**: When sampling a parameter (e.g., `super_schw_Amap`), remove it from `extra_args`. It can only appear in one place.

### 8.5 Build Artifacts After Code Changes

Always do a clean rebuild:

```bash
make clean
rm -rf classy classy_bak
make -j4
cd python && python3 setup.py build_ext --inplace && cd ..
```

Don't skip `make clean` — the `classy` binary and `python/classy.cpython-*.so` can cache stale code.

---

## 9. Appendix: Diagnostic Results

### 9.1 θ_s-Matched Fixed-Parameter Predictions

| Model | H₀_phys | H₀_local | σ₈ | S₈ | S₈_local |
|-------|---------|----------|------|------|----------|
| ΛCDM | 67.2 | 67.2 | 0.808 | 0.828 | 0.828 |
| Model A (unsq) | 67.2 | 67.2 | 0.752 | 0.801 | 0.801 |
| Model B (unsq, ODE dilution) | 67.2 | ~71.1 | 0.752 | 0.801 | ~0.745 |
| Model B (unsq, ODE no-decay) | 67.2 | ~71.6 | 0.752 | 0.801 | ~0.739 |
| Model B (unsq, accum) | 67.2 | ~72.6 | 0.752 | 0.801 | ~0.720 |

### 9.2 Observational Targets

- **SH0ES**: H₀ = 73.04 ± 1.04 km/s/Mpc
- **Planck CMB**: H₀ = 67.4 ± 0.5 km/s/Mpc
- **DES Y3**: S₈ = 0.776 ± 0.017
- **Planck S₈**: 0.834 ± 0.016

### 9.3 MCMC Performance Summary (unsquared, ieff_type=4)

| Run | Δχ² vs ΛCDM | H₀_local | S₈ | SH0ES χ² |
|-----|-------------|----------|------|----------|
| Model B A2 (ODE dilution) | −47 | 72.0 | 0.779 | 1.0 |
| Model B A2 (ODE no-decay) | −49 | 73.0 | 0.780 | 0.0 |
| Model B A-free (ODE dilution) | −47 | 73.1 | 0.779 | 0.0 |
| Model A | −31 | 68.9 | 0.777 | 15.0 |
| ΛCDM A2 (ODE, no interaction) | — | 72.3 | 0.800 | 0.3 |

### 9.4 Key Physical Results

1. **CPL surrogate**: Background indistinguishable from ΛCDM (|δH/H| < 0.01%). Entire S₈ shift is from perturbation interaction.
2. **f_clust**: Data independently finds f_clust ≈ 0.03 (consistent with zero). f_clust = 1 catastrophically worsens fit by +146 χ².
3. **I_eff exponent**: Squared is the minimum working exponent (fourth power too weak). Unsquared (linear) is preferred — stronger S₈ suppression, better fits, simpler theory.
4. **ODE dilution vs no-decay**: The three X₀ modes (accumulator, no-decay, full ODE) give X₀ = 0.035, 0.027, 0.022 respectively. σ₈ and S₈ are identical — only the H₀ mapping differs. The no-decay ODE has a clean physical argument: source weighted by horizon growth rate, but accumulated entropy doesn't dilute because thermodynamic production is irreversible. The full ODE is the most conservative (standard perturbation theory with no additional assumptions). MCMC will discriminate.
5. **β stability**: At unsquared coupling, β = 1/12 is already optimal. β-free runs show minimal movement from the theory prediction.
