#!/usr/bin/env python3
# frsd_transfer_v3.py -- publication figure: direct-computation curves
# (v2 engine) with the original figure's survey compilation & styling.
# Run from ~/class_holo_test:  python3 frsd_transfer_v3.py
import numpy as np
import matplotlib.pyplot as plt
from classy import Class

base = {
    'omega_b': 0.02237, 'omega_cdm': 0.12, 'h': 0.68,
    'A_s': 2.0989032e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
    'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
    'output': 'mPk,dTk,vTk', 'P_k_max_1/Mpc': 10.0, 'z_max_pk': 3.5,
    'gauge': 'newtonian',
}
holo = dict(base); holo.update({'interaction_beta': 1.0/12.0,
    'interaction_ieff_type': 4, 'f_clust': 0.0})
ZP = np.linspace(0.05, 2.2, 30)
KH = [0.02, 0.05, 0.10, 0.20]

DATA = [  # (z, fs8, err_lo, err_hi, survey, marker) -- verbatim compilation
    (0.067, 0.423, 0.055, 0.055, "6dFGS", "o"),
    (0.15,  0.490, 0.145, 0.145, "SDSS MGS", "s"),
    (0.38,  0.497, 0.045, 0.045, "BOSS DR12", "D"),
    (0.51,  0.458, 0.038, 0.038, "BOSS DR12", "D"),
    (0.61,  0.436, 0.034, 0.034, "BOSS DR12", "D"),
    (0.70,  0.473, 0.041, 0.041, "eBOSS LRG", "^"),
    (0.85,  0.315, 0.095, 0.095, "eBOSS ELG", "v"),
    (1.48,  0.462, 0.045, 0.045, "eBOSS QSO", "p"),
    (0.60,  0.550, 0.120, 0.120, "VIPERS", "*"),
    (0.86,  0.400, 0.110, 0.110, "VIPERS", "*"),
    (1.36,  0.482, 0.116, 0.116, "FastSound", "h"),
]
COL = {'6dFGS': '#4daf4a', 'SDSS MGS': '#984ea3', 'BOSS DR12': '#ff7f00',
       'eBOSS LRG': '#a65628', 'eBOSS ELG': '#f781bf', 'eBOSS QSO': '#999999',
       'VIPERS': '#377eb8', 'FastSound': '#e41a1c'}

def rho_interp(c, key):
    bg = c.get_background(); zb = bg['z']
    return lambda z: np.interp(z, zb[::-1], bg[key][::-1])

def run(params, label):
    c = Class(); c.set(params); c.compute()
    h = c.h(); rc, rb = rho_interp(c, '(.)rho_cdm'), rho_interp(c, '(.)rho_b')
    frs8, fds8 = [], []
    for z in ZP:
        s8 = c.sigma(8.0/h, z); Hc = c.Hubble(z)/(1+z)
        wc, wb = rc(z), rb(z)
        fr, fd = [], []
        tr = c.get_transfer(z); kk = tr['k (h/Mpc)'] * h
        trp = c.get_transfer(z+0.01); trm = c.get_transfer(z-0.01)
        wcp, wbp = rc(z+0.01), rb(z+0.01); wcm, wbm = rc(z-0.01), rb(z-0.01)
        kkp = trp['k (h/Mpc)'] * h; kkm = trm['k (h/Mpc)'] * h
        for kh in KH:
            k = kh*h
            g  = lambda t, K, n: np.interp(k, K, t[n])
            dm  = (wc*g(tr, kk, 'd_cdm') + wb*g(tr, kk, 'd_b'))/(wc+wb)
            tm  = (wc*g(tr, kk, 't_cdm') + wb*g(tr, kk, 't_b'))/(wc+wb)
            dmp = (wcp*g(trp, kkp, 'd_cdm') + wbp*g(trp, kkp, 'd_b'))/(wcp+wbp)
            dmm = (wcm*g(trm, kkm, 'd_cdm') + wbm*g(trm, kkm, 'd_b'))/(wcm+wbm)
            fr.append(-tm/(Hc*dm))
            fd.append(-(1+z)*(np.log(abs(dmp))-np.log(abs(dmm)))/0.02)
        frs8.append(np.mean(fr)*s8); fds8.append(np.mean(fd)*s8)
    c.struct_cleanup(); c.empty()
    print(f"{label}: done")
    return np.array(frs8), np.array(fds8)

print("Computing (30 redshifts x 2 cosmologies, a few minutes)...")
fL, dL = run(base, "LCDM"); fH, dH = run(holo, "Holographic")

fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6.6, 8.0),
    gridspec_kw={"height_ratios": [3, 1, 1]}, layout="constrained")
for a in ax:
    a.axvspan(0.02, 0.63, alpha=0.18, color='#e8c26e',
              label=r'over-capacity window ($\mathcal{S}>1$)' if a is ax[0] else None)
ax[0].plot(ZP, fL, '-',  color='#2166ac', lw=2, label=r'$\Lambda$CDM $f_{\rm RSD}\sigma_8$')
ax[0].plot(ZP, fH, '--', color='firebrick', lw=2, label=r'Model B $f_{\rm RSD}\sigma_8$ ($\beta=1/12$)')
ax[0].plot(ZP, dH, ':',  color='firebrick', lw=2, label=r'Model B $f_\delta\sigma_8$')
seen = set()
for z, y, el, eu, sv, mk in DATA:
    ax[0].errorbar(z, y, yerr=[[el], [eu]], fmt=mk, color=COL[sv], ms=6,
                   capsize=2, lw=1, label=(sv if sv not in seen else None))
    seen.add(sv)
ax[0].set_ylabel(r'$f\sigma_8(z)$'); ax[0].set_ylim(0.22, 0.65)
ax[0].legend(ncol=2, fontsize=7.5, frameon=True, loc='upper right')
ax[1].plot(ZP, fH/fL, '--', color='firebrick', lw=2)
ax[1].axhline(1, color='gray', lw=0.6, ls=':')
ax[1].set_ylim(0.972, 1.006)
ax[1].set_ylabel(r'Model B / $\Lambda$CDM')
ax[1].annotate(r'$1.9\%$ suppression at $z=0.1$', xy=(0.95, 0.982),
               fontsize=8, color='firebrick')
ax[2].plot(ZP, dH/fH, '-', color='k', lw=1.5)
ax[2].axhline(1, color='gray', lw=0.6, ls=':')
ax[2].set_ylim(0.70, 1.03); ax[2].set_ylabel(r'$f_\delta/f_{\rm RSD}$')
ax[2].annotate(r'density--velocity split: $26\%$ at $z\simeq0.1$',
               xy=(0.75, 0.78), fontsize=8)
ax[2].set_xlabel(r'Redshift $z$'); ax[2].set_xlim(0.0, 2.2)
fig.savefig('fsigma8_v3.pdf'); print("saved fsigma8_v3.pdf")
