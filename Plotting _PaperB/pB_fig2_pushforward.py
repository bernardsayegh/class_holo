#!/usr/bin/env python3
# Paper B Figure 2: joint (H0_phys, X_H) + prior-blind push-forward.
# Reads freeze_samples.npz (written by xh_posterior_freeze.py). Instant.
#   python3 pB_fig2_pushforward.py [freeze_samples.npz]
import sys, numpy as np, matplotlib.pyplot as plt
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "freeze_samples.npz")
XH, H0P = d["XH"], d["H0P"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
a1.scatter(H0P, XH, s=4, alpha=0.4)
a1.set_xlabel(r'$H_{0,\rm phys}$ [km/s/Mpc]'); a1.set_ylabel(r'$X_H$')
Hlad = H0P*np.exp(2*XH)
a2.hist(H0P, bins=30, alpha=0.5, density=True, label=r'$H_{0,\rm phys}$')
a2.hist(Hlad, bins=30, alpha=0.5, density=True,
        label=r'push-forward $H_{0,\rm phys}e^{2X_H}$')
a2.set_xlabel(r'$H_0$ [km/s/Mpc]'); a2.legend(frameon=False)
fig.savefig("pB_fig2_pushforward.pdf"); print("saved pB_fig2_pushforward.pdf")
