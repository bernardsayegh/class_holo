#!/usr/bin/env python3
# frsd_transfer_v2.py -- per-k validation + Figure 12 revision.
# Adds over v1: f_delta(k) directly from density transfers (no sigma8
# proxy), coded k/k_eq filter factor F(k), per-k residual table
# R(k,z) = f_RSD - f_delta(k) - F(k)*Q/(H rho), and the two/three-panel
# figure per the reviewer specification. Run from ~/class_holo_test.
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
ZT = [0.10, 0.30, 0.50, 0.70, 1.00, 1.50, 2.00]        # table redshifts
ZP = list(np.round(np.linspace(0.05, 2.0, 14), 3))      # plot redshifts
KH = [0.02, 0.05, 0.10, 0.20]
K_EQ = 0.073 * (base['omega_cdm'] + base['omega_b'])    # 1/Mpc
F = lambda k: (k/K_EQ)**2 / (1.0 + (k/K_EQ)**2)

def rho_interp(c, key):
    bg = c.get_background(); zb = bg['z']
    return lambda z: np.interp(z, zb[::-1], bg[key][::-1])

def species(c, z, k, wc, wb):
    tr = c.get_transfer(z); kk = tr['k (h/Mpc)'] * c.h()
    g = lambda n: np.interp(k, kk, tr[n])
    dm = (wc*g('d_cdm') + wb*g('d_b')) / (wc + wb)
    tm = (wc*g('t_cdm') + wb*g('t_b')) / (wc + wb)
    return dm, tm

def run(params, label, zs):
    c = Class(); c.set(params); c.compute()
    h = c.h(); rc, rb = rho_interp(c, '(.)rho_cdm'), rho_interp(c, '(.)rho_b')
    print(f"\n{label}  (sigma8(0) = {c.sigma8():.4f})")
    out = {}
    for z in zs:
        s8 = c.sigma(8.0/h, z); Hc = c.Hubble(z)/(1+z)
        wc, wb = rc(z), rb(z)
        Om = c.Om_m(z); QHr = (4.5*(1/12.0)*(1-0.75*Om)*(1-Om)
                               if 'interaction_beta' in params else 0.0)
        dz = 0.01
        rows = []
        for kh in KH:
            k = kh*h
            dm, tm = species(c, z, k, wc, wb)
            dp, _ = species(c, z+dz, k, rc(z+dz), rb(z+dz))
            dmn, _ = species(c, z-dz, k, rc(z-dz), rb(z-dz))
            fdk = -(1+z)*(np.log(abs(dp)) - np.log(abs(dmn)))/(2*dz)
            frsd = -tm/(Hc*dm)
            R = frsd - fdk - F(k)*QHr
            rows.append((kh, fdk, frsd, F(k), R))
        out[z] = (s8, QHr, rows)
    c.struct_cleanup(); c.empty()
    return out

print("="*78); print("PER-K VALIDATION TABLE"); print("="*78)
L = run(base, "LCDM", ZT); H = run(holo, "Holographic (beta=1/12)", ZT)
print(f"\nk_eq = {K_EQ:.5f}/Mpc   R(k,z) = f_RSD - f_delta(k) - F(k)*Q/Hrho")
print("  z    k[h/Mpc]   f_delta(k)  f_RSD(k)    F(k)      R (holo)     R (LCDM)")
for z in ZT:
    for (kh, fd, fr, Fk, R), (_, _, _, _, RL) in zip(H[z][2], L[z][2]):
        print(f" {z:4.2f}   {kh:5.2f}     {fd:8.4f}   {fr:8.4f}   {Fk:6.4f}   {R:+.6f}    {RL:+.6f}")

print("\nBuilding figure (14 redshifts x 2 cosmologies)...")
Lp = run(base, "LCDM (plot grid)", ZP); Hp = run(holo, "Holo (plot grid)", ZP)
def curves(D, zs):
    frs8, fds8 = [], []
    for z in zs:
        s8, _, rows = D[z]
        frs8.append(np.mean([r[2] for r in rows]) * s8)
        fds8.append(np.mean([r[1] for r in rows]) * s8)
    return np.array(frs8), np.array(fds8)
fL, dL = curves(Lp, ZP); fH, dH = curves(Hp, ZP)
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6.5, 8.2),
    gridspec_kw={"height_ratios": [3, 1, 1]}, layout="constrained")
ax[0].plot(ZP, fL, 'b-',  lw=2, label=r'$\Lambda$CDM $f_{\rm RSD}\sigma_8$')
ax[0].plot(ZP, fH, 'r--', lw=2, label=r'Holo $f_{\rm RSD}\sigma_8$')
ax[0].plot(ZP, dH, 'r:',  lw=2, label=r'Holo $f_\delta\sigma_8$')
ax[0].set_ylabel(r'$f\sigma_8$'); ax[0].legend(frameon=False)
ax[1].plot(ZP, fH/fL, 'k-', lw=1.5); ax[1].axhline(1, color='gray', lw=0.6)
ax[1].set_ylim(0.975, 1.005)
ax[1].set_ylabel(r'$\frac{(f_{\rm RSD}\sigma_8)_{\rm holo}}{(f_{\rm RSD}\sigma_8)_{\Lambda}}$')
ax[2].plot(ZP, dH/fH, 'k-', lw=1.5); ax[2].axhline(1, color='gray', lw=0.6)
ax[2].set_ylabel(r'$f_\delta/f_{\rm RSD}$'); ax[2].set_xlabel(r'$z$')
fig.savefig('fsigma8_v2.pdf'); print("saved fsigma8_v2.pdf")
