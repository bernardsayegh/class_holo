#!/usr/bin/env python3
"""
cosmo_stats_all_models.py
Compute cosmological statistics for ΛCDM and all holographic models (A, B, C, D).

This script:
1. Solves for h values that match the Planck θs constraint
2. Computes σ8, S8, Ωm, H0, H0_local for each model
3. Compares results in a summary table

Models:
- ΛCDM:    β=0, no interaction
- Model A: β=1/12, interaction only (no mapping, no reservoir)
- Model B: β=1/12, interaction + mapping (Amap=2, no reservoir)
- Model C: β=1/12, interaction + mapping + reservoir (Amap=2, amp=1)
- Model D: β=1/12, interaction + reservoir only (amp=1, no mapping)

Usage:
    python3 cosmo_stats_all_models.py
"""
import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np

# -----------------------------
# Targets / controls
# -----------------------------
THETA_S_TARGET = 1.040423   # target is 100*theta_s (Planck 2018)
H_BRACKET = (0.55, 0.80)
BETA = 1.0/12.0             # Holographic coupling
PKMAX = 10.0
# -----------------------------

# Base cosmological parameters
BASE_PARAMS = {
    "output": "mPk",
    "P_k_max_1/Mpc": PKMAX,
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

# Model configurations
MODELS = {
    "ΛCDM": {
        # No holographic parameters
    },
    "Model A": {
        "interaction_beta": BETA,
        "f_clust": 0.0,
        # SCR off (or on but inert)
    },
    "Model B": {
        "interaction_beta": BETA,
        "f_clust": 0.0,
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 0.0,       # No reservoir
        "super_schw_Amap": 2.0,      # Mapping on
        "super_schw_deltaS": 0.03,
        "super_schw_gamma": 2.0,
        "super_schw_no_mapping": 0,
    },
    "Model C": {
        "interaction_beta": BETA,
        "f_clust": 0.0,
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 1.0,       # Reservoir active
        "super_schw_Amap": 2.0,      # Mapping on
        "super_schw_deltaS": 0.03,
        "super_schw_gamma": 2.0,
        "super_schw_no_mapping": 0,
    },
    "Model D": {
        "interaction_beta": BETA,
        "f_clust": 0.0,
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 1.0,       # Reservoir active
        "super_schw_Amap": 2.0,
        "super_schw_deltaS": 0.03,
        "super_schw_gamma": 2.0,
        "super_schw_no_mapping": 1,  # Mapping OFF
    },
}


def compute_stats(h: float, model_params: dict):
    """Compute cosmological statistics for given h and model parameters."""
    p = dict(BASE_PARAMS)
    p["h"] = float(h)
    p.update(model_params)

    c = Class()
    c.set(p)
    c.compute()

    # Get basic derived parameters
    theta = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
    H0 = 100.0 * c.h()
    Om = float(c.Omega_m())
    sig8 = float(c.sigma8())
    S8 = sig8 * np.sqrt(Om / 0.3)

    # Try to get H0_local (only available if SCR is enabled)
    try:
        H0_local = c.get_current_derived_parameters(["H0_local"])["H0_local"]
    except:
        H0_local = H0  # Falls back to H0 if not available

    c.struct_cleanup()
    c.empty()
    return {
        "theta": float(theta),
        "H0": float(H0),
        "H0_local": float(H0_local),
        "Omega_m": Om,
        "sigma8": sig8,
        "S8": float(S8),
    }


def solve_h_for_theta(model_params: dict, target: float, bracket=(0.55, 0.80)):
    """Find h that gives target theta_s using Brent's method."""
    a, b = bracket

    def f(h):
        stats = compute_stats(h, model_params)
        return stats["theta"] - target

    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        raise RuntimeError(
            f"Root not bracketed. fa={fa:+.6e}, fb={fb:+.6e}. "
            "Try widening H_BRACKET."
        )

    return brentq(f, a, b, xtol=1e-8, rtol=1e-12, maxiter=200)


def run_model(name: str, model_params: dict):
    """Run a model: solve for h, compute stats, print results."""
    print(f"\n{'='*60}")
    print(f"{name}")
    print(f"{'='*60}")
    
    # Print model configuration
    if model_params:
        print("Configuration:")
        for k, v in model_params.items():
            print(f"  {k}: {v}")
    else:
        print("Configuration: Standard ΛCDM (no holographic parameters)")
    
    # Solve for h
    print(f"\nSolving for h to match θs = {THETA_S_TARGET}...")
    h_sol = solve_h_for_theta(model_params, THETA_S_TARGET, H_BRACKET)
    
    # Compute full stats
    stats = compute_stats(h_sol, model_params)
    
    print(f"\nResults:")
    print(f"  h (solved)    = {h_sol:.6f}")
    print(f"  100*theta_s   = {stats['theta']:.6f} (target: {THETA_S_TARGET})")
    print(f"  H0 (physical) = {stats['H0']:.2f} km/s/Mpc")
    print(f"  H0_local      = {stats['H0_local']:.2f} km/s/Mpc")
    print(f"  Omega_m       = {stats['Omega_m']:.5f}")
    print(f"  sigma8        = {stats['sigma8']:.5f}")
    print(f"  S8            = {stats['S8']:.5f}")
    
    return {"h": h_sol, **stats}


def print_summary_table(results: dict):
    """Print a comparison table of all models."""
    print("\n")
    print("=" * 90)
    print("SUMMARY TABLE")
    print("=" * 90)
    
    # Header
    print(f"{'Model':<12} {'h':>8} {'H0':>8} {'H0_loc':>8} {'Ωm':>8} {'σ8':>8} {'S8':>8}")
    print("-" * 90)
    
    for name, stats in results.items():
        print(f"{name:<12} {stats['h']:>8.5f} {stats['H0']:>8.2f} {stats['H0_local']:>8.2f} "
              f"{stats['Omega_m']:>8.5f} {stats['sigma8']:>8.5f} {stats['S8']:>8.5f}")
    
    print("-" * 90)
    
    # Comparison with targets
    print("\nComparison with observational targets:")
    print(f"  SH0ES H0:      73.04 ± 1.04 km/s/Mpc")
    print(f"  Planck H0:     67.4 ± 0.5 km/s/Mpc")
    print(f"  DES Y3 S8:     0.776 ± 0.017")
    print(f"  Planck S8:     0.834 ± 0.016")
    
    # Check which model resolves tensions
    print("\nTension resolution:")
    for name, stats in results.items():
        h0_ok = abs(stats['H0_local'] - 73.04) < 2.0
        s8_ok = abs(stats['S8'] - 0.776) < 0.03
        status = []
        if h0_ok:
            status.append("H0 ✓")
        if s8_ok:
            status.append("S8 ✓")
        if not status:
            status.append("—")
        print(f"  {name:<12}: {', '.join(status)}")


if __name__ == "__main__":
    print("=" * 90)
    print("Holographic Dark Energy: All Models Comparison")
    print("=" * 90)
    print(f"Target: 100*theta_s = {THETA_S_TARGET}")
    print(f"Coupling: β = 1/12 = {BETA:.8f}")
    print(f"h bracket: [{H_BRACKET[0]}, {H_BRACKET[1]}]")
    
    results = {}
    
    for name, params in MODELS.items():
        try:
            results[name] = run_model(name, params)
        except Exception as e:
            print(f"\n*** ERROR running {name}: {e}")
            continue
    
    if results:
        print_summary_table(results)
    
    print("\nDone.")
