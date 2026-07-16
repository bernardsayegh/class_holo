#!/usr/bin/env python3
# Post-process the frozen non-ladder chain: per-sample X_H, z_entry, z_exit,
# H0_phys; report <X_H>, sigma(X_H), Corr(X_H, H0_phys); emit Paper-B
# Figures 1 & 2.   Run from ~/class_holo_test after the freeze chain:
#   python3 xh_posterior_freeze.py chains_freeze/modelB_noshoes.1.txt
import sys, numpy as np, matplotlib.pyplot as plt
from classy import Class

root = sys.argv[1]          # root or any .N.txt of the chain
if root.endswith(".txt"):
    root = root.rsplit(".", 2)[0]
import glob
files = sorted(glob.glob(root + ".[0-9].txt"))
assert files, f"no chain files match {root}.N.txt"
names = open(files[0]).readline().lstrip("#").split()
col = {n: i for i, n in enumerate(names)}
parts = []
for f in files:
    d = np.loadtxt(f)
    parts.append(d[int(0.30*len(d)):])      # 30% burn-in per chain
raw = np.vstack(parts)
print(f"{len(files)} chains, {len(raw)} post-burn-in samples; columns: {names[:12]}...")
THIN = max(1, len(raw) // 400)
rows = raw[::THIN]
print(f"{len(rows)} thinned samples (thin={THIN})")

def one(pars):
    c = Class()
    c.set({'omega_b': pars['omega_b'], 'omega_cdm': pars['omega_cdm'],
           'H0': pars['H0'], 'tau_reio': pars['tau_reio'],
           'n_s': pars['n_s'], 'A_s': pars['A_s'],
           'YHe': 0.2454, 'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
           'interaction_beta': 0.0833333, 'interaction_ieff_type': 4,
           'f_clust': 0.0,
           'super_schwarzschild_correction': 'yes',
           'super_schw_amp': 0.0, 'super_schw_Amap': 2.0,
           'super_schw_deltaS': 0.01, 'super_schw_gamma': 2.0,
           'super_schw_ode': 0,
           'output': ''})   # background-only rerun; Pk settings irrelevant to X_H
    c.compute()
    bg = c.get_background(); z = bg['z']; a = 1/(1+z)
    OL = bg['(.)rho_lambda']/bg['(.)rho_crit']
    S = 4.5*OL*(1-OL)
    D = np.clip(1-1/np.maximum(S, 1e-12), 0, None)
    lna = np.log(a); o = np.argsort(lna)
    XH = np.trapz(D[o], lna[o])
    Xcum = np.concatenate([[0], np.cumsum(0.5*(D[o][1:]+D[o][:-1])*np.diff(lna[o]))])
    # S=1 crossings
    f = S[o]-1; zs = z[o]; idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    cross = [zs[i]-f[i]*(zs[i+1]-zs[i])/(f[i+1]-f[i]) for i in idx]
    H0p = c.Hubble(0)*299792.458   # computed physical H(0), NOT input h
    c.struct_cleanup(); c.empty()
    return XH, (max(cross) if cross else np.nan), (min(cross) if cross else np.nan), H0p, (zs, S[o], D[o], Xcum)

recs, curves = [], []
for r in rows:
    p = {k: r[col[k]] for k in ('omega_b','omega_cdm','H0','tau_reio','n_s')}
    p['A_s'] = r[col['A_s']] if 'A_s' in col else 1e-10*np.exp(r[col['logA']])
    XH, zin, zout, H0p, cur = one(p)
    recs.append((XH, zin, zout, H0p)); curves.append(cur)
R = np.array(recs); XH, ZIN, ZOUT, H0P = R.T
# ── self-certification against the production pipeline's own mapping ──
if 'H0_local' in col:
    H0L = rows[:, col['H0_local']]
    dev = np.abs(H0P*np.exp(2*XH) - H0L)/H0L
    print(f"CROSS-CHECK vs stored H0_local: max rel. dev = {dev.max():.2e}"
          f"  (mean {dev.mean():.2e}) -- should be ~1e-3 or below")
print(f"<X_H> = {XH.mean():.5f}   sigma(X_H) = {XH.std():.5f}")
print(f"z_entry = {ZIN.mean():.3f} +/- {ZIN.std():.3f}   z_exit = {ZOUT.mean():.3f} +/- {ZOUT.std():.3f}")
print(f"Corr(X_H, H0_phys) = {np.corrcoef(XH, H0P)[0,1]:+.3f}")
np.savez("freeze_samples.npz", XH=XH, H0P=H0P, ZIN=ZIN, ZOUT=ZOUT)

# ── Figure 1: S(z), Delta(z), cumulative X_H(z) with posterior band ──
zg = np.linspace(0, 3, 300)
def band(ax, idx, lab):
    def sinterp(c):
        zi = np.asarray(c[0]); yi = np.asarray(c[idx]); o = np.argsort(zi)
        return np.interp(zg, zi[o], yi[o])
    ys = np.array([sinterp(c) for c in curves])
    lo, md, hi = np.percentile(ys, [16, 50, 84], axis=0)
    ax.fill_between(zg, lo, hi, alpha=0.3); ax.plot(zg, md, lw=1.5)
    ax.set_ylabel(lab)
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6, 7), layout="constrained")
band(ax[0], 1, r'$\mathcal{S}(z)$'); ax[0].axhline(1, color='gray', ls=':', lw=0.7)
band(ax[1], 2, r'$\Delta(z)$')
band(ax[2], 3, r'$X_H(<z)$'); ax[2].set_xlabel(r'$z$'); ax[2].invert_xaxis()
fig.savefig("pB_fig1_capacity_history.pdf"); print("saved pB_fig1")

# ── Figure 2: joint (H0_phys, X_H) + prior-blind push-forward ──
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
a1.scatter(H0P, XH, s=4, alpha=0.4)
a1.set_xlabel(r'$H_{0,\rm phys}$'); a1.set_ylabel(r'$X_H$')
Hlad = H0P*np.exp(2*XH)
a2.hist(H0P, bins=30, alpha=0.5, density=True, label=r'$H_{0,\rm phys}$')
a2.hist(Hlad, bins=30, alpha=0.5, density=True,
        label=r'push-forward $H_{0,\rm phys}e^{2X_H}$')
a2.set_xlabel(r'$H_0$ [km/s/Mpc]'); a2.legend(frameon=False)
fig.savefig("pB_fig2_pushforward.pdf"); print("saved pB_fig2")
