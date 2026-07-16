#!/usr/bin/env python3
# Paper-B Figures 3 & 4 (data-free):  python3 fig_pB_diagrams.py
import numpy as np, matplotlib.pyplot as plt
X0 = 0.0359; u = 5*X0/np.log(10)         # 0.07796 mag
# ── Fig 3: native plane (Delta M_B, 5 Delta a_B) ──
fig, ax = plt.subplots(figsize=(5.4, 5.4), layout="constrained")
L = 0.18; t = np.linspace(-L, L, 2)
ax.plot(t, t,  'firebrick', lw=2, label=r'common gate: $\Delta M_B = 5\Delta a_B$')
ax.plot(t, -t, 'gray', lw=1.5, ls='--',
        label=r'photometric null: $\Delta M_B = -5\Delta a_B$')
ax.scatter([u], [u], s=70, color='firebrick', zorder=5)
ax.annotate(rf'$(\,{u:.5f},\,{u:.5f}\,)$ at $X_H={X0}$',
            xy=(u, u), xytext=(u-0.155, u+0.02), fontsize=9)
ax.annotate('', xy=(0.115, -0.115), xytext=(0.05, -0.05),
            arrowprops=dict(arrowstyle='->', color='gray'))
ax.text(0.075, -0.108, r'$D_{\rm com}=\Delta M_B-5\Delta a_B$', fontsize=9, color='gray')
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel(r'$5\,\Delta a_B$  [mag]'); ax.set_ylabel(r'$\Delta M_B$  [mag]')
ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect('equal')
ax.legend(loc='lower right', fontsize=8, frameon=False)
fig.savefig("pB_fig3_native_plane.pdf"); print("saved pB_fig3")
