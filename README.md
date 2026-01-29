# Holographic Dark Energy Implementation in CLASS

Modified version of [CLASS](https://github.com/lesgourg/class_public) implementing holographic dark energy with geometric dark sector interactions.

> **Holographic Equilibrium: Simultaneous Resolution of the H₀ and S₈ Tensions**

## Quick Start

```bash
# Compile
make clean && make class

# Run Model B (resolves both S₈ and H₀)
./class explanatory_holo_B.ini
```

## Key Parameters

| Parameter | Model A | Model B | Description |
|-----------|---------|---------|-------------|
| `interaction_beta` | 0.0833 | 0.0833 | Coupling strength (β=1/12) |
| `f_clust` | 0.0 | 0.0 | Clustering fraction |
| `super_schwarzschild_correction` | no | yes | Enable SCR subsystem |
| `super_schw_Amap` | - | 2.0 | H₀ mapping amplitude |
| `super_schw_amp` | - | 0.0 | Reservoir strength |

## Results

- **S₈ tension resolved**: σ₈ = 0.749 ± 0.005 (vs ΛCDM: 0.798)
- **H₀ tension resolved**: H₀_local = 73.30 km/s/Mpc (Model B)
- **Fit improvement**: Δχ² = -55.0 vs ΛCDM

## Branch

Use `option-c-scr-q` branch for the holographic implementation.

## Citation

[Paper citation to be added]
