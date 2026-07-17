#!/usr/bin/env python3
"""Paper B: screening response of the explicit two-leg scheduler.

The default output name is kept unchanged:
``paper_B_screening_scheduler.pdf``.

The plotted curves are

    x_inf(f) = f(1+2f)/3,
    x_N(f)   = x_inf(f) - 2f(1-f)/[3(N-1)],

and the same-slot activity observable ``f``.  The finite-population curve is
an exact ordered-pair count for fixed active population in the continuum
plotting approximation.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def screening_large_n(f: np.ndarray) -> np.ndarray:
    return f * (1.0 + 2.0 * f) / 3.0


def screening_fixed_population(f: np.ndarray, n_slots: int) -> np.ndarray:
    if n_slots < 2:
        raise ValueError("N must be at least 2")
    return screening_large_n(f) - 2.0 * f * (1.0 - f) / (3.0 * (n_slots - 1))


def make_figure(*, n_slots: int = 20, points: int = 500, output: str | Path = "paper_B_screening_scheduler.pdf") -> Path:
    if points < 20:
        raise ValueError("points must be at least 20")
    f = np.linspace(0.0, 1.0, points)
    x_inf = screening_large_n(f)
    x_n = screening_fixed_population(f, n_slots)

    # Exact endpoint checks keep the plotting script tied to the equations.
    if not np.allclose([x_inf[0], x_inf[-1]], [0.0, 1.0], atol=1e-14):
        raise RuntimeError("large-N screening endpoint check failed")
    if np.any(x_n < -1e-13) or np.any(x_n > 1.0 + 1e-13):
        raise RuntimeError("finite-population screening left the physical interval")

    fig, ax = plt.subplots(figsize=(5.6, 4.2), layout="constrained")
    ax.plot(f, x_inf, linewidth=2.0,
            label=r"two-leg large-$N$: $x_\infty=f(1+2f)/3$")
    ax.plot(f, x_n, linewidth=1.5, linestyle="--",
            label=rf"fixed population, $N={n_slots}$ (exact)")
    ax.plot(f, f, linewidth=1.2, linestyle=":",
            label=r"same-slot activity $f$")
    ax.set_xlabel(r"adjoint activity fraction $f$")
    ax.set_ylabel(r"normalized response $x$")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.03)
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", dest="n_slots", type=int, default=20,
                        help="slot count used for the exact fixed-population curve (default: 20)")
    parser.add_argument("--points", type=int, default=500)
    parser.add_argument("--output", type=Path, default=Path("paper_B_screening_scheduler.pdf"))
    args = parser.parse_args()
    out = make_figure(n_slots=args.n_slots, points=args.points, output=args.output)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
