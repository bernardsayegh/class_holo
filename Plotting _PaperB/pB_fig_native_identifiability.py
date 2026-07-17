#!/usr/bin/env python3
"""Paper B: native identifiability geometry (illustrative, not a fit).

The ordinary ladder fixes a ridge of constant ``x_T+x_R``.  The common-gate
condition ``x_T=x_R`` selects one point only after an additional native
reference or closure condition is supplied.  The output filename is kept
unchanged: ``paper_B_native_identifiability.pdf``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_figure(*, summed_loading: float = 0.0718,
                output: str | Path = "paper_B_native_identifiability.pdf") -> Path:
    if not np.isfinite(summed_loading):
        raise ValueError("summed_loading must be finite")

    margin = max(0.035, 0.42 * abs(summed_loading))
    lo = min(-0.025, -margin)
    hi = max(0.095, summed_loading + margin)
    x_t = np.linspace(lo, hi, 500)
    ridge = summed_loading - x_t
    common = x_t
    x_star = summed_loading / 2.0

    fig, ax = plt.subplots(figsize=(5.3, 5.0), layout="constrained")
    ax.plot(x_t, ridge, linewidth=2.3,
            label=rf"illustrative ladder ridge: $x_T+x_R={summed_loading:.4f}$")
    ax.plot(x_t, common, linewidth=1.8, linestyle="--",
            label=r"common-gate condition: $x_T=x_R$")
    ax.scatter([x_star], [x_star], s=62, zorder=5)
    ax.annotate("intersection after an\nadditional native condition",
                xy=(x_star, x_star), xytext=(x_star + 0.014, x_star + 0.030),
                fontsize=8.5, ha="left",
                arrowprops=dict(arrowstyle="->", linewidth=0.8))

    # The ladder-null motion in the (x_T,x_R) plane preserves the sum.
    p0 = np.array([x_star - 0.018, x_star + 0.018])
    p1 = np.array([x_star + 0.018, x_star - 0.018])
    ax.annotate("", xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle="<->", linewidth=1.0))
    ax.text(x_star + 0.020, x_star - 0.024,
            r"null direction $(+1,-1)$", fontsize=8.2, rotation=-38)

    ax.axhline(0.0, linewidth=0.7)
    ax.axvline(0.0, linewidth=0.7)
    ax.set_xlabel(r"temporal loading $x_T$")
    ax.set_ylabel(r"ruler loading $x_R$")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.20)
    ax.legend(frameon=False, fontsize=8.2, loc="lower left")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sum", dest="summed_loading", type=float, default=0.0718,
                        help="illustrative value of x_T+x_R (default: 0.0718)")
    parser.add_argument("--output", type=Path, default=Path("paper_B_native_identifiability.pdf"))
    args = parser.parse_args()
    out = make_figure(summed_loading=args.summed_loading, output=args.output)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
