# class_holo

Modified [CLASS](https://github.com/lesgourg/class_public) Boltzmann solver implementing holographic dark energy with negative-feedback braking at the apparent horizon.

Companion code for:  
**"Holographic Braking in ΛCDM: Resolving the H₀ and S₈ Tensions via Apparent-Horizon Feedback"**  
B. Sayegh (2026)

---

## What this does

The standard CLASS solver is extended with an interaction between dark energy and the apparent horizon, governed by a single coupling constant β = 1/12 derived from Schwarzschild–de Sitter stability. The modification:

- Suppresses late-time clustering (σ₈ ↓) via holographic drag on dark matter perturbations
- Generates a local Hubble constant H₀,loc > H₀,CMB through geometric mapping at the super-Schwarzschild boundary
- Preserves w = −1 (no phantom crossing, no extra dynamical degrees of freedom)

Model B (β = 1/12, A_map = 2) simultaneously resolves the H₀ tension (H₀,loc = 73.30 km/s/Mpc) and S₈ tension (S₈ = 0.771) with Δχ² = −55.0 relative to ΛCDM.

---

## Repository structure

```
class_holo/
├── source/                       # Modified CLASS C source (background + perturbations)
├── python/                       # Python wrapper (classy)
├── configs/                      # Cobaya YAML configuration files
│   ├── lcdm_shoes.input.yaml
│   ├── lcdm_no_h0prior.input.yaml
│   ├── lcdm_trgb.input.yaml
│   ├── modelA_fixed.input.yaml
│   ├── modelA_f1.input.yaml
│   ├── modelA_f_free.input.yaml
│   ├── modelA_beta_free.input.yaml
│   ├── modelB_Amap2.input.yaml
│   ├── modelB_Amap1.input.yaml
│   ├── modelB_Amap_free.input.yaml
│   ├── modelB_beta_free.input.yaml
│   ├── modelB_no_h0prior.input.yaml
│   ├── modelB_trgb.input.yaml
│   ├── modelC_theory.input.yaml
│   ├── modelC_amp_free.input.yaml
│   ├── modelD_theory.input.yaml
│   ├── modelD_amp_free.input.yaml
│   ├── acdm_Amap_free.input.yaml
│   └── acdm_mapping_Amap2.input.yaml
├── likelihoods/                  # Custom likelihood modules
│   ├── desy3_likelihood.py
│   ├── desy3_likelihood_h0local.py
│   └── shoes_h0local.py
├── scripts/                      # Plotting and analysis scripts
│   ├── plot_evolution_combined.py
│   ├── plot_holographic_data_bracket_expanded.py
│   ├── plot_spectra_combined.py
│   ├── plot_contours_combined.py
│   ├── plot_posteriors.py
│   ├── plot_S8_combined.py
│   ├── plot_triangle_S8_Om_H0.py
│   ├── plot_fsigma8_comparison.py
│   ├── plot_desi_bao_combined.py
│   └── cosmo_stats_thetas_matched.py   # Validation diagnostic
└── README.md
```

---

## Installation

### Build CLASS

```bash
git clone https://github.com/bernardsayegh/class_holo.git
cd class_holo
make clean
rm -rf classy classy_bak
make -j4
cd python && python3 setup.py build_ext --inplace && cd ..
```

### Verify

Quick smoke test:

```bash
python3 -c "
import sys; sys.path.insert(0, 'python')
from classy import Class
c = Class()
c.set({'output':'mPk','P_k_max_1/Mpc':10,'omega_b':0.02242,'omega_cdm':0.11933,'h':0.6766,'n_s':0.9665,'A_s':2.1e-9,'tau_reio':0.054})
c.compute()
print(f'OK: sigma8={c.sigma8():.4f}')
"
```

### Validation

A more thorough check that holographic braking and geometric mapping are working correctly. Fixes ω_b, ω_cdm to Planck 2018 values and solves for h such that 100×θₛ = 1.040423 in each model. (Note: h values differ from MCMC best-fits because the MCMC varies all parameters simultaneously.)

```bash
python3 scripts/cosmo_stats_thetas_matched.py
```

Expected output:

```
ΛCDM (β=0)
  h (solved)   = 0.67164851
  100*theta_s  = 1.040423
  H0           = 67.1649 km/s/Mpc
  H0_local     = 67.1649 km/s/Mpc
  Omega_m      = 0.31565160
  sigma8       = 0.80761068
  S8           = 0.82841018

Model A (β=1/12)
  h (solved)   = 0.66783350
  100*theta_s  = 1.040423
  H0           = 66.7833 km/s/Mpc
  H0_local     = 67.2795 km/s/Mpc
  Omega_m      = 0.32927198
  sigma8       = 0.74254936
  S8           = 0.77793281

Model B (β=1/12, Amap=2)
  h (solved)   = 0.66202759
  100*theta_s  = 1.040423
  H0           = 66.2028 km/s/Mpc
  H0_local     = 73.1651 km/s/Mpc
  Omega_m      = 0.31984210
  sigma8       = 0.73750851
  S8           = 0.76150756
```

What to check:

- **σ₈ drops ~8%** from ΛCDM → Model A/B (holographic drag working)
- **H₀,loc ≈ 73.2** for Model B (geometric mapping working)
- **H₀,loc ≈ H₀** for ΛCDM (no mapping, as expected)
- **100×θₛ = 1.040423** for all three (CMB acoustic scale preserved)
- **Model B has lower h but higher H₀,loc than Model A** — the mapping modifies the background expansion, so the solver compensates with lower h to preserve θₛ while the locally-measured value jumps up

### Python environment

```bash
python3 -m venv venv
source venv/bin/activate
pip install cobaya getdist matplotlib numpy scipy pyccl
```

### Install Planck and BAO data

```bash
cobaya-install cosmo --packages-path packages
```

### DES-Y3 data

The DES-Y3 2pt data file (`2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits`) is too large for GitHub. Download from the [DES Data Release](https://des.ncsa.illinois.edu/releases/y3a2) and place in the `likelihoods/` directory.

---

## Reproducing results

Each YAML file is self-contained. To reproduce any chain:

```bash
source venv/bin/activate
cobaya-run configs/<file>.input.yaml
```

### Configuration → Result mapping

| Configuration file | Produces |
|---|---|
| `lcdm_shoes.input.yaml` | ΛCDM baseline (Tables 1, 2) |
| `modelA_fixed.input.yaml` | Model A, β = 1/12 fixed (Table 2) |
| `modelA_f1.input.yaml` | Model A, f_clust = 1 variant (Table 2) |
| `modelA_f_free.input.yaml` | Model A, f_clust free (Table 2) |
| `modelA_beta_free.input.yaml` | Model A, β free (Table 2) |
| `modelB_Amap2.input.yaml` | Model B headline (Tables 1, 2, 3) |
| `modelB_Amap1.input.yaml` | Model B, A_map = 1 variant (Table 2) |
| `modelB_Amap_free.input.yaml` | Model B, A_map free (Table 2) |
| `modelB_beta_free.input.yaml` | Model B, β free (Table 2) |
| `modelC_theory.input.yaml` | Model C, theory prior (Table 2) |
| `modelC_amp_free.input.yaml` | Model C, amplitude free (Table 2) |
| `modelD_theory.input.yaml` | Model D, theory prior (Table 2) |
| `modelD_amp_free.input.yaml` | Model D, amplitude free (Table 2) |
| `acdm_Amap_free.input.yaml` | ΛCDM + mapping controls (Table 4) |
| `acdm_mapping_Amap2.input.yaml` | ΛCDM + mapping, A_map = 2 (Table 4) |
| `modelB_no_h0prior.input.yaml` | H₀ prior robustness (Table 5) |
| `modelB_trgb.input.yaml` | H₀ prior robustness (Table 5) |
| `lcdm_no_h0prior.input.yaml` | ΛCDM no-prior control (Table 5) |
| `lcdm_trgb.input.yaml` | ΛCDM TRGB control (Table 5) |

---

## Custom likelihoods

### `desy3_likelihood.py`
DES-Y3 cosmic shear (ξ±) using pyccl with NLA intrinsic alignments. Computes shear correlation functions via the BBKS transfer function; the clustering amplitude is set by σ₈ from modified CLASS. Switching to Eisenstein–Hu changes the DES χ² by ~2, well within sampling noise.

### `desy3_likelihood_h0local.py`
Local-frame variant for ΛCDM control runs. Inherits from the base DES likelihood but evaluates S₈ using the optically-mapped H₀,loc and corresponding Ωₘ,loc.

### `shoes_h0local.py`
Gaussian likelihood on H₀,loc (the locally-measured Hubble constant output by modified CLASS), centred on the SH0ES measurement.

---

## Plotting scripts

All scripts are standalone. Run from the repository root with CLASS on the Python path:

```bash
PYTHONPATH=python python3 scripts/<script>.py
```

Scripts that load MCMC chains (getdist) expect a `chains_schw/` directory with converged chains.

| Script | Figure | Description |
|---|---|---|
| `plot_evolution_combined.py` | Fig. 2 (top) | Density parameter evolution Ωₘ(z), ΩΛ(z) |
| `plot_holographic_data_bracket_expanded.py` | Fig. 2 (bottom-left) | Holographic data bracket with ΩΛ = 1/3, 2/3 roots |
| `plot_spectra_combined.py` | Fig. 3 | CMB TT/EE/lensing power spectra comparison |
| `plot_contours_combined.py` | Fig. 4 | 2D posterior contours (key parameter pairs) |
| `plot_posteriors.py` | Fig. 5 | 1D marginalised posteriors |
| `plot_S8_combined.py` | Fig. 6 | S₈ constraints across models and datasets |
| `plot_triangle_S8_Om_H0.py` | Fig. 7 | Corner plot: H₀, σ₈, Ωₘ, S₈ with tension bands |
| `plot_fsigma8_comparison.py` | Fig. 8 | fσ₈(z) comparison with RSD survey data |
| `plot_desi_bao_combined.py` | Fig. 9 | DESI Y1 BAO comparison (D_M/r_d, D_H/r_d) |

Fig. 2 (bottom-right) is rendered as inline TikZ in the paper source — no external script needed.

---

## Key CLASS parameters

The holographic modification is controlled by these parameters in the `.ini` or YAML files:

| Parameter | Description | Model B value |
|---|---|---|
| `interaction_beta` | Coupling constant β | 0.0833333 (= 1/12) |
| `f_clust` | Dark energy clustering fraction | 0.0 |
| `super_schwarzschild_correction` | Enable geometric mapping | `yes` |
| `super_schw_Amap` | Mapping amplitude | 2.0 |
| `super_schw_deltaS` | Transition width | 0.03 |
| `super_schw_gamma` | Transition steepness | 2.0 |

Setting `interaction_beta: 0` and `super_schwarzschild_correction: no` recovers standard ΛCDM.

---

## Citation

If you use this code, please cite:

```bibtex
@article{Sayegh2026,
    author  = {Sayegh, Bernard},
    title   = {Holographic Braking in {$\Lambda$CDM}: Resolving the {$H_0$} 
               and {$S_8$} Tensions via Apparent-Horizon Feedback},
    year    = {2026},
    note    = {Submitted to Physical Review D}
}
```

---

## License

The CLASS modifications are released under the same licence as [CLASS](https://github.com/lesgourg/class_public) (GNU GPLv3). Custom likelihoods and scripts are MIT-licensed.
