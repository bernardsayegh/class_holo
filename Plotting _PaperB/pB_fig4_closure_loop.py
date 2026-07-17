#!/usr/bin/env python3
"""Paper B Figure 4: bright-siren versus SN-ladder closure loop.

The two routes are drawn through their actual intermediate standards rather
than as arcs that bypass the GW and SN nodes.  The default output name is
unchanged.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np


def _arrow(ax: plt.Axes, p0: tuple[float, float], p1: tuple[float, float], *,
           connectionstyle: str = "arc3,rad=0.0", linewidth: float = 1.9,
           linestyle: str = "-") -> None:
    patch = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=12,
                            linewidth=linewidth, linestyle=linestyle, connectionstyle=connectionstyle)
    ax.add_patch(patch)


def make_figure(*, x_h: float = 0.0359,
                output: str | Path = "pB_fig4_closure_loop.pdf") -> Path:
    if not np.isfinite(x_h):
        raise ValueError("x_h must be finite")
    ratio = np.exp(2.0 * x_h)
    percent = 100.0 * (ratio - 1.0)
    mag = 10.0 * x_h / np.log(10.0)

    fig, ax = plt.subplots(figsize=(6.5, 4.75), layout="constrained")
    ax.axis("off")
    nodes = {
        "E": (0.08, 0.50),
        "GW": (0.47, 0.78),
        "SN": (0.47, 0.22),
        "O": (0.92, 0.50),
    }
    labels = {
        "E": "shared source /\nhost volume",
        "GW": "bright-siren\ndistance",
        "SN": "SN absolute\ncalibration",
        "O": "common\nobserver",
    }
    for key, (x, y) in nodes.items():
        ax.scatter([x], [y], s=900, facecolor="white", edgecolor="black", zorder=3)
        ax.text(x, y, labels[key], ha="center", va="center", fontsize=8.2, zorder=4)

    # Neutral gravitational-wave path: E -> GW -> O.
    _arrow(ax, nodes["E"], nodes["GW"], connectionstyle="arc3,rad=0.05")
    _arrow(ax, nodes["GW"], nodes["O"], connectionstyle="arc3,rad=-0.05")
    ax.text(0.49, 0.93, r"neutral route: $\Delta\mathbf{W}=(0,0)$",
            ha="center", fontsize=9.0)

    # Full common-gate electromagnetic ladder path: E -> SN -> O.
    _arrow(ax, nodes["E"], nodes["SN"], connectionstyle="arc3,rad=-0.05", linestyle="--")
    _arrow(ax, nodes["SN"], nodes["O"], connectionstyle="arc3,rad=0.05", linestyle="--")
    ax.text(0.49, 0.055, r"ladder route: $\Delta\mathbf{W}=(1,1)$",
            ha="center", fontsize=9.0)

    ax.text(0.51, 0.50,
            rf"$\mathcal{{H}}_{{\rm GW-SN}}=2X_H={2*x_h:.4f}$" + "\n"
            + rf"route ratio $={ratio:.5f}$ ({percent:.2f}\%)" + "\n"
            + rf"$={mag:.5f}$ mag-equivalent",
            ha="center", va="center", fontsize=9.7,
            bbox=dict(boxstyle="round,pad=0.42", facecolor="white", edgecolor="0.5"))

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xh", type=float, default=0.0359)
    parser.add_argument("--output", type=Path, default=Path("pB_fig4_closure_loop.pdf"))
    args = parser.parse_args()
    out = make_figure(x_h=args.xh, output=args.output)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
