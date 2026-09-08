# consistency-audit branch

Diagnostic branch for the background/perturbation consistency questions raised in the
September 2026 review of Papers A and B. **No production result is changed by this branch:
every switch defaults to legacy behaviour**, and the legacy row of the audit matrix reproduces
the certified numbers (X_H = 0.03588, H0_phys = 68.487, S8 = 0.7745 at the fixed point).

## The problem, in one paragraph

The production code keeps rho_Lambda constant and injects Q into CDM (`background.c`: the Lambda
block sets `rho_lambda = Omega0_lambda*H0^2`; the derivative block has
`dy[rho_cdm] = -3 rho_cdm + Q_over_H`). The papers' displayed equations describe an
energy-conserving vacuum donor (`rho_Lambda' = -Q`). Two consequences follow for the code as
written: (1) with no compensating pressure, nabla_mu T^{mu nu} = 0 is violated, and the stored
`H' = -(3/2)(rho+p)a` is not the derivative of the integrated H(a) (6.4% error at low z);
(2) the accumulator gate uses the *reference* parabola S = (9/2)u(1-u), whereas the sweep
derivation evaluated on the actual background gives S_geom = 3 Omega_L (1+q_geom), which
nearly halves the exposure. The certified 0.03588 is therefore a hybrid quantity.

## Switches added (all default 0)

| input parameter | effect |
|---|---|
| `interaction_creation_pressure = 1` | adds Pi = -Q/(3H) to `p_tot` (Route B background completion). Restores the Bianchi identity; H(a) unchanged. `p_tot_prime` not updated (feeds only the unused PPF branch). |
| `interaction_gate_geometric = 1` | accumulator uses S_geom = 3 Omega_L^full (1+q_geom) with q_geom from the actual dlnH/dlna. |
| `interaction_cdm_closure = 1` | CDM continuity uses the creation-pressure closure (dPi=0, theta_c=0): delta' = -(1-Gamma/3)(theta+h'/2) - aH Gamma delta, unfiltered, no f_clust. Newtonian theta equation left legacy. |
| `interaction_vacuum_donor = 1` | Route A: rho_Lambda integrated with d rho_L/dlna = -Q/H. |

The kernel Q(rho_de, rho_b, rho_cdm) is factored into `holo_kernel_Q_over_H()` so the same
expression serves the derivative block, the creation pressure and the geometric gate.

## Audit matrix at the certified fixed point (`audit/results/audit_matrix_results.json`)

| config | H0_phys | X_H | H0_ladder | H' residual | sigma8 | S8 | max dTT (l<=50) |
|---|---|---|---|---|---|---|---|
| legacy | 68.487 | 0.03588 | 73.583 | 3.2e-2 | 0.7444 | 0.7745 | -- |
| creation_pressure | 68.487 | 0.03588 | 73.583 | 3.7e-4 | 0.7444 | 0.7745 | 5.1% (l=2) |
| geometric_gate | 68.487 | **0.01968** | **71.236** | 3.2e-2 | 0.7444 | 0.7745 | ~0 |
| cdm_closure | 68.487 | 0.03588 | 73.583 | 3.2e-2 | **0.7346** | **0.7643** | 9.1% (l=2) |
| closure+creation | 68.487 | 0.03588 | 73.583 | 3.7e-4 | 0.7346 | 0.7643 | 4.6% (l=2) |
| vacuum_donor | **39.56** | 0.000 | 39.56 | 1.5e-4 | 0.997 | 1.73 | 31% |

Reading:
- **Creation pressure** fixes the background consistency (H' residual down 85x) at zero cost to
  H(a), sigma8, S8. Its CMB signature is confined to l < 30 (late ISW); rough TT-only
  dchi2 against cosmic variance = 0.05. A joint-likelihood evaluation is still required.
- **Geometric gate** reproduces the reviewer's independent result (0.019624 / 71.22). Under the
  literal sweep derivation the ladder register lands ~2 sigma below R25. Keeping 0.03588 requires
  stating the reference parabola as a constitutive response law.
- **Closure** shifts sigma8 by -1.3% (the sub-horizon proxy predicted -1.9%). S8 = 0.7643 is
  still 0.66 sigma from DES Y3, but the posterior must be re-derived; old chains cannot be assumed.
- **Vacuum donor with this kernel is not viable**: Q ~ 4.5 beta I_eff rho_Lambda at early times
  (I_eff -> 1/4), so d rho_L/dlna ~ -0.094 rho_L over ~32 e-folds and Lambda drains by a factor 20.
  The reviewer's reference calculation started at a = 0.01 and saw only the last 4.6 e-folds.
  Route A needs a kernel that vanishes faster at early times before it can be tested.

## Running

```
make -j                                  # C library
cd python && python3 setup.py build_ext --inplace && cd ..   # wrapper, in place
python3 audit/experiments/run_audit_matrix.py      # ~3 min: the table above
python3 audit/experiments/tt_lowell_profile.py     # l-profile of the CMB shifts
python3 audit/diagnostics/geometric_sweep.py       # pure-python cross-check, no CLASS
```

Chains: any switch can be passed through cobaya's `extra_args` exactly like `interaction_beta`.

## What this branch does not do

It does not decide the physics. Route B (constant Lambda + creation pressure) needs a
justified perturbative closure -- the one implemented here is the simplest, not the only one --
and re-derived chains. Route A needs a modified kernel before it has a background at all. The
Paper B dust-clock action assumes the donor form (V_T = -Q) and must be re-derived under
whichever route is chosen. The papers remain paused until one route is chosen and re-fit.
