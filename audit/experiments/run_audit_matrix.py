#!/usr/bin/env python3
"""Consistency-audit matrix at the certified Paper-A fixed point.

Runs the same cosmology under each audit switch and reports, per configuration:
  H0_in, H0_phys, X_H (code accumulator), H0_ladder = H0_phys*exp(2 X_H),
  the H'-consistency residual (numerical dlnH/dlna vs the stored -(3/2)(rho+p)/rho),
  sigma8 and S8, and the max fractional change in C_l^TT (l<=50) vs legacy.

Run from a built class_holo checkout:   python3 audit/experiments/run_audit_matrix.py
"""
import sys, json, numpy as np
sys.path.insert(0, "python")
from classy import Class

BASE = {'omega_b': 0.022539, 'omega_cdm': 0.117409, 'H0': 67.6255, 'A_s': 2.096955e-9,
        'n_s': 0.970111, 'tau_reio': 0.056507, 'YHe': 0.2454, 'N_ur': 2.0328, 'N_ncdm': 1,
        'm_ncdm': 0.06, 'interaction_beta': 0.0833333, 'interaction_ieff_type': 4, 'f_clust': 0.0,
        # accumulator on with zero reservoir amplitude: records X_H and H0_local, no backreaction
        'super_schwarzschild_correction': 'yes', 'super_schw_amp': 0.0, 'super_schw_Amap': 2.0,
        'super_schw_deltaS': 0.01, 'super_schw_gamma': 2.0, 'super_schw_ode': 0,
        'output': 'tCl,lCl,mPk', 'lensing': 'yes', 'l_max_scalars': 2500, 'P_k_max_1/Mpc': 3.0}
C = 299792.458

CONFIGS = [
    ("legacy",            {}),
    ("creation_pressure", {'interaction_creation_pressure': 1}),
    ("geometric_gate",    {'interaction_gate_geometric': 1}),
    ("geom+creation",     {'interaction_gate_geometric': 1, 'interaction_creation_pressure': 1}),
    ("cdm_closure",       {'interaction_cdm_closure': 1}),
    ("closure+creation",  {'interaction_cdm_closure': 1, 'interaction_creation_pressure': 1}),
    ("vacuum_donor",      {'interaction_vacuum_donor': 1}),
]

def run(extra):
    c = Class(); c.set({**BASE, **extra}); c.compute()
    bg = c.get_background()
    z = np.asarray(bg['z']); H = np.asarray(bg['H [1/Mpc]'])
    rho = np.asarray(bg['(.)rho_tot']); p = np.asarray(bg['(.)p_tot'])
    lna = -np.log1p(z); o = np.argsort(lna)
    dlnH_num = np.gradient(np.log(H[o]), lna[o])
    dlnH_stored = -1.5 * (rho[o] + p[o]) / rho[o]
    late = z[o] < 1.0
    hres = float(np.max(np.abs(dlnH_num[late] - dlnH_stored[late])))
    XH = float(np.asarray(bg['X_schw'])[o][-1]) if 'X_schw' in bg else float('nan')
    H0p = float(c.Hubble(0)) * C
    out = dict(H0_in=BASE['H0'], H0_phys=H0p, XH=XH, H0_ladder=H0p * np.exp(2 * XH),
               Hprime_resid=hres, sigma8=float(c.sigma8()), Omega_m=float(c.Omega_m()))
    out['S8'] = out['sigma8'] * np.sqrt(out['Omega_m'] / 0.3)
    cl = c.lensed_cl(2500); out['_tt'] = np.asarray(cl['tt'])
    rl = np.asarray(bg['(.)rho_lambda'])[o]; out['rhoL_ratio_today_vs_early'] = float(rl[-1] / rl[0])
    c.struct_cleanup(); c.empty()
    return out

res = {}
for name, extra in CONFIGS:
    try:
        res[name] = run(extra); print(f"ran {name}")
    except Exception as e:
        print(f"FAILED {name}: {str(e)[:200]}"); res[name] = None

ref = res.get("legacy")
print(f"\n{'config':<18}{'H0_phys':>9}{'X_H':>10}{'H0_ladder':>10}{'H'' resid':>11}{'sigma8':>9}{'S8':>8}{'dTT(l<=50)':>12}{'rhoL(0)/rhoL(i)':>16}")
print("-" * 105)
for name, _ in CONFIGS:
    r = res[name]
    if r is None: print(f"{name:<18}  FAILED"); continue
    dtt = float(np.max(np.abs(r['_tt'][2:51] / ref['_tt'][2:51] - 1))) if ref else float('nan')
    print(f"{name:<18}{r['H0_phys']:>9.3f}{r['XH']:>10.5f}{r['H0_ladder']:>10.3f}{r['Hprime_resid']:>11.2e}"
          f"{r['sigma8']:>9.4f}{r['S8']:>8.4f}{dtt:>12.2e}{r['rhoL_ratio_today_vs_early']:>16.4f}")
json.dump({k: ({kk: vv for kk, vv in v.items() if not kk.startswith('_')} if v else None) for k, v in res.items()},
          open("audit_matrix_results.json", "w"), indent=2)
print("\nwrote audit_matrix_results.json")
print("\nREFERENCE VALUES: legacy X_H 0.03588 (certified), H0_phys 68.487 at this point; geometric gate ~0.0196 (reviewer 0.019624);")
print("closure sub-horizon estimate sigma8 ratio ~0.981; donor H0_phys/H0_in ~0.876 with rho_L depleted.")
