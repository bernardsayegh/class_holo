#!/usr/bin/env python3
"""Paper B Figure 3: native common-gate and photometric-null directions.

The script is analytic and data-free.  Its filename and default output name
are unchanged.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_figure(*, x_h: float = 0.0359,
                output: str | Path = "pB_fig3_native_plane.pdf") -> Path:
    if not np.isfinite(x_h):
        raise ValueError("x_h must be finite")
    kappa_x = 5.0 * x_h / np.log(10.0)
    lim = max(0.18, 2.15 * abs(kappa_x))
    t = np.linspace(-lim, lim, 400)

    fig, ax = plt.subplots(figsize=(5.35, 5.15), layout="constrained")
    ax.plot(t, t, linewidth=2.1,
            label=r"common gate: $\Delta M_B=5\Delta a_B$")
    ax.plot(t, -t, linewidth=1.5, linestyle="--",
            label=r"photometric null: $\Delta M_B=-5\Delta a_B$")
    ax.scatter([kappa_x], [kappa_x], s=62, zorder=5)
    ax.annotate(rf"$X_H={x_h:.4f}$" + "\n" + rf"$({kappa_x:.5f},{kappa_x:.5f})$ mag",
                xy=(kappa_x, kappa_x), xytext=(-0.155, kappa_x + 0.035),
                fontsize=8.5,
                arrowprops=dict(arrowstyle="->", linewidth=0.8))

    # Show the orthogonal displacement measured by D_com.
    q0 = np.array([0.055, 0.055])
    q1 = np.array([0.100, 0.010])
    ax.annotate("", xy=q1, xytext=q0,
                arrowprops=dict(arrowstyle="->", linewidth=1.0))
    ax.text(0.082, 0.018,
            r"off-gate: $D_{\rm com}=\Delta M_B-5\Delta a_B$",
            fontsize=7.9, rotation=-36, ha="center")

    ax.axhline(0.0, linewidth=0.6)
    ax.axvline(0.0, linewidth=0.6)
    ax.set_xlabel(r"$5\,\Delta a_B$ [mag]")
    ax.set_ylabel(r"$\Delta M_B$ [mag]")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.18)
    ax.legend(loc="lower right", fontsize=8.0, frameon=False)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xh", type=float, default=0.0359)
    parser.add_argument("--output", type=Path, default=Path("pB_fig3_native_plane.pdf"))
    args = parser.parse_args()
    out = make_figure(x_h=args.xh, output=args.output)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
