#!/usr/bin/env python3
# Paper B Figure 1: S(z), Delta(z), cumulative X_H(<z) with posterior band.
# Standalone. Run from ~/class_holo_test:
#   python3 pB_fig1_capacity_history.py chains_test/modelB_A2_noprior_v5_acc.1.txt
# Uses N_CURVES thinned samples (band is rigidity-thin; 60 suffices, ~2 min).
import sys, glob, numpy as np, matplotlib.pyplot as plt
from classy import Class

N_CURVES = 60
root = sys.argv[1]
if root.endswith(".txt"):
    root = root.rsplit(".", 2)[0]
files = sorted(glob.glob(root + ".[0-9].txt"))
names = open(files[0]).readline().lstrip("#").split()
col = {n: i for i, n in enumerate(names)}
parts = [np.loadtxt(f)[int(0.30*len(np.loadtxt(f))):] for f in files]
raw = np.vstack(parts)
rows = raw[::max(1, len(raw)//N_CURVES)]
print(f"{len(rows)} samples for curves")

def curves(r):
    c = Class()
    c.set({'omega_b': r[col['omega_b']], 'omega_cdm': r[col['omega_cdm']],
           'H0': r[col['H0']], 'tau_reio': r[col['tau_reio']],
           'n_s': r[col['n_s']],
           'A_s': r[col['A_s']] if 'A_s' in col else 1e-10*np.exp(r[col['logA']]),
           'YHe': 0.2454, 'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
           'interaction_beta': 0.0833333, 'interaction_ieff_type': 4,
           'f_clust': 0.0, 'super_schwarzschild_correction': 'yes',
           'super_schw_amp': 0.0, 'super_schw_Amap': 2.0,
           'super_schw_deltaS': 0.01, 'super_schw_gamma': 2.0,
           'super_schw_ode': 0, 'output': ''})
    c.compute()
    bg = c.get_background(); z = bg['z']; a = 1/(1+z)
    OL = bg['(.)rho_lambda']/bg['(.)rho_crit']
    S = 4.5*OL*(1-OL)
    D = np.clip(1-1/np.maximum(S, 1e-12), 0, None)
    lna = np.log(a); o = np.argsort(lna)
    Xc = np.concatenate([[0], np.cumsum(0.5*(D[o][1:]+D[o][:-1])*np.diff(lna[o]))])
    c.struct_cleanup(); c.empty()
    return z[o], S[o], D[o], Xc

zg = np.linspace(0, 3, 300)
all_c = [curves(r) for r in rows]
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6, 7), layout="constrained")
for idx, lab in [(1, r'$\mathcal{S}(z)$'), (2, r'$\Delta(z)$'), (3, r'$X_H(<z)$')]:
    ys = np.array([np.interp(zg, c[0], c[idx]) for c in all_c])
    lo, md, hi = np.percentile(ys, [16, 50, 84], axis=0)
    a_ = ax[idx-1]
    a_.fill_between(zg, lo, hi, alpha=0.3); a_.plot(zg, md, lw=1.5)
    a_.set_ylabel(lab)
ax[0].axhline(1, color='gray', ls=':', lw=0.7)
ax[2].set_xlabel(r'$z$'); ax[2].invert_xaxis()
fig.savefig("pB_fig1_capacity_history.pdf"); print("saved pB_fig1_capacity_history.pdf")
