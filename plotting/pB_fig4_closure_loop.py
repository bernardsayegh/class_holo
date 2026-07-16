#!/usr/bin/env python3
# Paper-B Figures 3 & 4 (data-free):  python3 fig_pB_diagrams.py
import numpy as np, matplotlib.pyplot as plt
X0 = 0.0359; u = 5*X0/np.log(10)         # 0.07796 mag
# ── Fig 4: bright-siren vs SN closure loop ──
fig, ax = plt.subplots(figsize=(6.2, 4.6), layout="constrained")
ax.axis('off')
P = {'E': (0.08, 0.5), 'GW': (0.5, 0.82), 'SN': (0.5, 0.18), 'O': (0.92, 0.5)}
for k, (x, y) in P.items():
    ax.scatter([x], [y], s=900, facecolor='white', edgecolor='k', zorder=3)
    ax.text(x, y, {'E':'source\nevent','GW':'siren path\n(GW amplitude)',
                   'SN':'ladder path\n(SN calibration)','O':'observer'}[k],
            ha='center', va='center', fontsize=8, zorder=4)
def arrow(a, b, dy, lab, col):
    (x1,y1),(x2,y2) = P[a], P[b]
    ax.annotate('', xy=(x2, y2+dy), xytext=(x1, y1+dy),
                arrowprops=dict(arrowstyle='->', lw=2, color=col,
                                connectionstyle=f"arc3,rad={0.25 if dy>0 else -0.25}"))
    ax.text((x1+x2)/2, y1+dy+(0.14 if dy>0 else -0.14), lab,
            ha='center', fontsize=9, color=col)
arrow('E', 'O',  0.32, r'$\Delta\mathbf{W}=(0,0)$: no standardization crossing', 'steelblue')
arrow('E', 'O', -0.32, r'$\Delta\mathbf{W}=(1,1)$: two-endpoint gate engaged', 'firebrick')
ax.text(0.5, 0.5, r'$\mathcal{H}_{\rm GW-SN}=2X_H=7.44\%$'
        '\n' r'$=0.15591$ mag-equivalent', ha='center', fontsize=11,
        bbox=dict(boxstyle='round', fc='ivory', ec='gray'))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
fig.savefig("pB_fig4_closure_loop.pdf"); print("saved pB_fig4")
