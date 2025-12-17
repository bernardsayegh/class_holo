# Holographic Dark Drag: Modified CLASS Code

This repository contains a modified version of the CLASS Boltzmann code implementing the holographic dark energy interaction model described in [your paper reference].

## Quick Start
```bash
# Clone the repository
git clone https://github.com/bernardsayegh/class_holo.git
cd class_holo

# Build CLASS
make clean
make -j4

# Install Python wrapper
cd python
pip install . --break-system-packages
cd ..

# Test the installation
python3 << 'PYEOF'
import sys
sys.path.insert(0, 'python')
from classy import Class
cosmo = Class()
cosmo.set({
    'output': 'mPk',
    'h': 0.6766,
    'omega_b': 0.02242,
    'omega_cdm': 0.11933,
    'Omega_Lambda': 0.0,
    'Omega_k': 0.0,
    'w0_fld': -0.962,
    'wa_fld': -0.07,
    'interaction_beta': 0.5,
    'interaction_area_dilution': 1,
    'interaction_use_ah_filter': 1,
    'f_clust': 0.0,
})
cosmo.compute()
print(f"Omega_m = {cosmo.Omega_m():.4f}")
print(f"sigma8 = {cosmo.sigma8():.4f}")
S8 = cosmo.sigma8() * (cosmo.Omega_m()/0.3)**0.5
print(f"S8 = {S8:.4f}")
cosmo.struct_cleanup()
cosmo.empty()
PYEOF
```

## New Parameters

The holographic interaction is controlled by the following parameters:

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| `interaction_beta` | Fundamental coupling strength β | 0.0 | 0.5 |
| `interaction_area_dilution` | Enable area dilution (0/1) | 0 | 1 |
| `interaction_use_ah_filter` | Enable horizon filter (0/1) | 0 | 1 |
| `f_clust` | Clustering fraction (0=homogeneous, 1=clustered) | 0.0 | 0.0 |

## Configuration for Paper Results

To reproduce the results in the paper, use:
```python
params = {
    'output': 'mPk',
    'P_k_max_1/Mpc': 10.0,
    'h': 0.6766,
    'omega_b': 0.02242,
    'omega_cdm': 0.11933,
    'n_s': 0.9665,
    'A_s': 2.1e-9,
    'Omega_Lambda': 0.0,      # IMPORTANT: Disable cosmological constant
    'Omega_k': 0.0,           # Flat universe
    'w0_fld': -0.962,         # Thawing dark energy
    'wa_fld': -0.07,
    'interaction_beta': 0.5,  # Holographic coupling
    'interaction_area_dilution': 1,
    'interaction_use_ah_filter': 1,
    'f_clust': 0.0,
}
```

## Expected Results

With the above configuration:

| Observable | ΛCDM | Holographic | 
|------------|------|-------------|
| σ₈ | 0.821 | 0.731 |
| S₈ | 0.834 | 0.775 |
| Ω_m | 0.310 | 0.337 |

## Modified Files

The holographic interaction is implemented in:

- `source/background.c` - Modified CDM continuity equation
- `source/perturbations.c` - Modified perturbation equations
- `include/background.h` - New parameter declarations
- `source/input.c` - Parameter parsing

## Physics Summary

The model modifies the CDM continuity equation:
```
dρ_cdm/dt = -3Hρ_cdm(1 - β_eff·Ω_de)
```

where:
- β_eff = β_fund / min(area_ratio, 4)
- area_ratio = (τ·a·H)²
- The interaction sources CDM from the cosmological horizon

This creates an attractor at Ω_m → 1/3, Ω_de → 2/3.

## Citation

If you use this code, please cite:
[Your paper reference]

## License

Based on CLASS (https://github.com/lesgourg/class_public)
