# class_holo

Modified [CLASS](https://github.com/lesgourg/class_public) Boltzmann solver implementing holographic dark energy–dark matter interaction via apparent-horizon thermodynamics.

Companion code for:  
**"Holographic Braking in ΛCDM: Horizon Kinematics and the Resolution of Cosmological Tensions"**  
B. Sayegh (2026)

---

## What this does

The standard CLASS solver is extended with a dark energy → dark matter interaction governed by a coupling constant β = 1/12 — the exact ratio of the Pendry-limited quantum channel throughput (c⁵/24G) to the classical geometric capacity (c⁵/2G) of the holographic screen. The modification:

- Suppresses late-time clustering (σ₈ ↓) via holographic drag on dark matter perturbations
- Generates a local Hubble constant H₀,loc > H₀,phys through a super-Schwarzschild Weyl-sector mapping
- Preserves w = −1 (no phantom crossing, no extra dynamical degrees of freedom)

Model B (β = 1/12, A_map = 2) simultaneously resolves the H₀ tension (H₀,loc = 72.65 ± 0.18 km/s/Mpc) and S₈ tension (S₈ = 0.778 ± 0.005) with Δχ² = −58.4 relative to ΛCDM using Planck 2018 + DESI DR2 + Pantheon+ + DES Y3 + SH0ES.

---

## Repository structure

```
class_holo/
├── README.md
├── IMPLEMENTATION.md
├── configs/                      # 15 Cobaya YAML configuration files
│   ├── lcdm.yaml
│   ├── modelB_A2.yaml
│   └── ...
├── cobaya/likelihoods/shoes_h0local/
│   ├── shoes_h0local.py
│   ├── desy3_likelihood.py
│   └── desy3_likelihood_h0local.py
├── plotting/                     # Figure scripts + validation scripts
├── source/                       # Modified CLASS C source
├── include/                      # Modified CLASS headers
├── python/                       # Python wrapper (classy)
├── scripts/                      # Original CLASS example scripts
└── ...                           # Standard CLASS directories
```

---

## Quick start

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
python3 plotting/cosmo_stats_thetas_matched.py
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

### ΛCDM+mapping control test

Confirms that the geometric mapping does nothing without holographic braking:

```bash
python3 plotting/check_acdm_mapping.py
```

Compares ΛCDM, ΛCDM+mapping (A_map=2, β=0), Model A (β=1/12), and Model B (A_map=2, β=1/12) at fixed h, printing H₀, H₀,loc, σ₈, S₈, S₈,loc, Ωₘ, Ωₘ,loc, and X₀. What to check:

- **ΛCDM+mapping** produces X₀ = 0 and H₀,loc = H₀ — mapping infrastructure present but inert without β
- **Model A** suppresses σ₈ but H₀,loc ≈ H₀ — braking without mapping
- **Model B** suppresses σ₈ *and* boosts H₀,loc — both mechanisms active
- **S₈,loc < S₈** for Model B because Ωₘ,loc < Ωₘ (same ωₘ, larger H₀,loc denominator)

---

## Python environment

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

The DES-Y3 2pt data file (`2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits`) is too large for GitHub. Download from the [DES Data Release](https://des.ncsa.illinois.edu/releases/y3a2) and place in the `cobaya/likelihoods/shoes_h0local/` directory.

---

## Reproducing results

Each YAML file is self-contained. To reproduce any chain:

```bash
source venv/bin/activate
cobaya-run configs/<file>.yaml
```

All chains converge to Gelman–Rubin R − 1 < 0.01 with the first 30% discarded as burn-in.

### Configuration → Result mapping

| Configuration file | Produces |
|---|---|
| `lcdm.yaml` | ΛCDM baseline |
| `lcdm_A2.yaml` | Frame control (β = 0, A_map = 2) |
| `lcdm_Afree.yaml` | ΛCDM + A_map free |
| `modelA.yaml` | Model A (β = 1/12) |
| `modelA_bfree.yaml` | Model A, β free |
| `modelA_fc1.yaml` | Model A, f_clust = 1 |
| `modelA_fcfree.yaml` | Model A, f_clust free |
| `modelB_A1.yaml` | Model B (A_map = 1) |
| `modelB_A2.yaml` | **Model B (A_map = 2), headline** |
| `modelB_A2_noprior.yaml` | Model B, no SH0ES prior |
| `modelB_A2_ode1.yaml` | Model B, dilution ODE variant |
| `modelB_A2_trgb.yaml` | Model B, TRGB calibration |
| `modelB_Afree.yaml` | Model B, A_map free |
| `modelC_Afree_Sfree.yaml` | Model C (A_map, S₀ free) |
| `modelD_Sfree.yaml` | Model D (S₀ free) |

### Likelihood stack

- `planck_2018_highl_plik.TTTEEE`
- `planck_2018_lowl.TT`
- `planck_2018_lowl.EE`
- `planck_2018_lensing.clik`
- `bao.desi_dr2` (custom; 7-bin BAO likelihood)
- `sn.pantheonplus`
- `desy3_likelihood.DESY3` (custom; physical S₈)
- `shoes_h0local.SH0ES_H0local` (custom; evaluates against H₀,loc)

---

## Custom likelihoods

### `desy3_likelihood.py`
DES-Y3 cosmic shear (ξ±) using pyccl with NLA intrinsic alignments. Computes shear correlation functions via the BBKS transfer function; the clustering amplitude is set by σ₈ from modified CLASS. Evaluates against physical S₈ (not the mapped S₈,loc).

### `desy3_likelihood_h0local.py`
Local-frame variant for ΛCDM control runs. Inherits from the base DES likelihood but evaluates S₈ using the optically-mapped H₀,loc and corresponding Ωₘ,loc.

### `shoes_h0local.py`
Gaussian likelihood on H₀,loc (the locally-measured Hubble constant output by modified CLASS), centred on the SH0ES measurement (H₀ = 73.04 ± 1.04 km/s/Mpc).

---

## Plotting scripts

All scripts are standalone. Run from the repository root with CLASS on the Python path:

```bash
PYTHONPATH=python python3 plotting/<script>.py
```

Scripts that load MCMC chains (getdist) expect a `chains_schw/` directory with converged chains.

| Script | Figure | Description |
|---|---|---|
| `plot_evolution_combined.py` | evolution | Density parameter evolution Ωₘ(z), ΩΛ(z) |
| `plot_holographic_data_bracket_expanded.py` | sweep_rate | Holographic data bracket with ΩΛ = 1/3, 2/3 roots |
| `plot_spectra_combined.py` | spectra | CMB TT/EE/lensing and matter P(k) comparison |
| `plot_contours_combined.py` | posteriors | 2D posterior contours (key parameter pairs) |
| `plot_posteriors.py` | posteriors | 1D marginalised posteriors |
| `plot_S8_combined.py` | S8_barchart | S₈ constraints across models and datasets |
| `plot_triangle_S8_Om_H0.py` | triangle | Corner plot: H₀, σ₈, Ωₘ, S₈ with tension bands |
| `plot_fsigma8_comparison.py` | fsigma8 | fσ₈(z) comparison with RSD survey data |
| `plot_desi_bao_combined.py` | desi_bao | DESI DR2 BAO comparison (D_M/r_d, D_H/r_d) |

The logic flow diagram (Fig. 1) is rendered as inline TikZ in the paper source — no external script needed.

---

## Key CLASS parameters

The holographic modification is controlled by these parameters in the `.ini` or YAML files:

| Parameter | Description | Model B value |
|---|---|---|
| `interaction_beta` | Coupling constant β | 0.0833333 (= 1/12) |
| `f_clust` | Clustering fraction of injected CDM | 0.0 |
| `super_schwarzschild_correction` | Enable geometric mapping | `yes` |
| `super_schw_Amap` | Mapping amplitude A_map | 2.0 |
| `super_schw_deltaS` | Transition smoothing width ΔS | 0.01 |
| `super_schw_amp` | Reservoir amplitude (Models C/D) | 0.0 |

Setting `interaction_beta: 0` and `super_schwarzschild_correction: no` recovers standard ΛCDM.

---

## Citation

If you use this code, please cite:

```bibtex
@article{Sayegh2026,
    author  = {Sayegh, Bernard},
    title   = {Holographic Braking in {$\Lambda$CDM}: Horizon Kinematics
               and the Resolution of Cosmological Tensions},
    year    = {2026},
    note    = {Submitted to Physical Review D}
}
```

---

## License

The CLASS modifications are released under the same licence as [CLASS](https://github.com/lesgourg/class_public) (GNU GPLv3). Custom likelihoods and scripts are MIT-licensed.
