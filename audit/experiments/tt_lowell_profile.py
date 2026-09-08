#!/usr/bin/env python3
"""Low-ell TT profile of each audit switch vs legacy, against per-multipole cosmic variance.
Run from the repo root after building the wrapper in place:  python3 audit/experiments/tt_lowell_profile.py"""
import sys, numpy as np; sys.path.insert(0, "python"); from classy import Class
exec(open("audit/experiments/run_audit_matrix.py").read().split("CONFIGS = [")[0])
B = {k: v for k, v in BASE.items() if k != 'P_k_max_1/Mpc'}; B['output'] = 'tCl,lCl'
def tt(extra):
    c = Class(); c.set({**B, **extra}); c.compute(); r = np.asarray(c.lensed_cl(2500)['tt']); c.struct_cleanup(); c.empty(); return r
ref = tt({}); l = np.arange(2, 2501); cv = np.sqrt(2 / (2 * l + 1))
for name, extra in [("creation_pressure", {'interaction_creation_pressure': 1}), ("cdm_closure", {'interaction_cdm_closure': 1}),
                    ("closure+creation", {'interaction_cdm_closure': 1, 'interaction_creation_pressure': 1})]:
    d = tt(extra)[2:] / ref[2:] - 1; i = np.argmax(abs(d[:49]))
    print(f"{name:<18} max|dTT|={abs(d[:49]).max():.3f} at l={l[i]}  l=30:{abs(d[28]):.4f}  l=100:{abs(d[98]):.4f}"
          f"  rough TT-only dchi2 vs cosmic variance = {np.sum((d / cv) ** 2):.2f}")
