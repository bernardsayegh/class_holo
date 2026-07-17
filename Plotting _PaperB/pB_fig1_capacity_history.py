#!/usr/bin/env python3
"""Paper B Figure 1: capacity history from a frozen non-ladder chain.

This is the standalone wrapper around the corrected reconstruction routines in
``xh_posterior_freeze.py``.  The script name and default output filename are
unchanged.

Example
-------

    python3 pB_fig1_capacity_history.py chains_freeze/modelB_noshoes.1.txt
"""
from __future__ import annotations

import argparse
from pathlib import Path

from xh_posterior_freeze import FreezeConfig, evaluate_chain, plot_capacity_history


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain", type=Path, help="chain root or any root.N.txt file")
    parser.add_argument("--burn-in", type=float, default=0.30)
    parser.add_argument("--curves", type=int, default=60,
                        help="number of posterior-resampled curves (default: 60)")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--beta", type=float, default=1.0 / 12.0)
    parser.add_argument("--delta-s", type=float, default=0.01)
    parser.add_argument("--hard-gate", action="store_true")
    parser.add_argument("--z-max", type=float, default=3.0)
    parser.add_argument("--output", type=Path, default=Path("pB_fig1_capacity_history.pdf"))
    args = parser.parse_args()

    config = FreezeConfig(
        burn_in=args.burn_in,
        n_samples=args.curves,
        seed=args.seed,
        beta=args.beta,
        delta_s=args.delta_s,
        hard_gate=args.hard_gate,
        z_plot_max=args.z_max,
    )
    _, _, _, _, records = evaluate_chain(args.chain, config)
    out = plot_capacity_history(records, output=args.output, z_max=args.z_max)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
