#!/usr/bin/env python3
"""Paper B Figure 2: sample-level prior-blind push-forward.

This version keeps the actual frozen posterior resamples as the evidence
layer.  The near-zero width of ``X_H`` is made visible by plotting the
mean-subtracted residual in automatically chosen scientific units rather than
replacing the samples by Gaussian summary curves.

The right panel likewise uses empirical sample histograms for
``H0_phys`` and ``H0_phys*exp(2*X_H)``.  Summary moments are annotated, but no
Gaussian shape is imposed.

The script name is unchanged.  By default it writes both
``pB_fig2_pushforward.pdf`` (legacy project name) and
``pB_fig2_pushforward_v9.pdf`` (the current Paper-B TeX name).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _finite_paired(x_h: np.ndarray, h_phys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_h = np.asarray(x_h, dtype=float).reshape(-1)
    h_phys = np.asarray(h_phys, dtype=float).reshape(-1)
    if x_h.size != h_phys.size:
        raise ValueError("XH and H0P must have the same number of samples")
    keep = np.isfinite(x_h) & np.isfinite(h_phys)
    x_h, h_phys = x_h[keep], h_phys[keep]
    if x_h.size < 3:
        raise ValueError("freeze file must contain at least three finite paired samples")
    return x_h, h_phys


def _quantiles(array: np.ndarray) -> tuple[float, float, float]:
    q16, q50, q84 = np.percentile(array, [16.0, 50.0, 84.0])
    return float(q16), float(q50), float(q84)


def _residual_scale(sigma: float) -> tuple[float, int]:
    """Choose 10**exponent so the residual standard deviation is O(1)."""
    if not np.isfinite(sigma) or sigma <= 0.0:
        return 1.0, 0
    exponent = max(0, int(np.floor(-np.log10(sigma))))
    return float(10.0**exponent), exponent


def samples_from_npz(path: str | Path) -> dict[str, float | np.ndarray]:
    with np.load(path) as data:
        if "XH" not in data or "H0P" not in data:
            raise KeyError("freeze file must contain arrays named XH and H0P")
        x_h, h_phys = _finite_paired(data["XH"], data["H0P"])
    h_ladder = h_phys * np.exp(2.0 * x_h)
    sigma_x = float(np.std(x_h, ddof=1))
    corr = float(np.corrcoef(x_h, h_phys)[0, 1]) if sigma_x > 0.0 else float("nan")
    return {
        "XH": x_h,
        "H0P": h_phys,
        "H0L": h_ladder,
        "mu_x": float(np.mean(x_h)),
        "sigma_x": sigma_x,
        "corr": corr,
        "mu_phys": float(np.mean(h_phys)),
        "sigma_phys": float(np.std(h_phys, ddof=1)),
        "mu_ladder": float(np.mean(h_ladder)),
        "sigma_ladder": float(np.std(h_ladder, ddof=1)),
    }


def _common_histogram_edges(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    combined = np.concatenate([a, b])
    edges = np.histogram_bin_edges(combined, bins="fd")
    if edges.size < 12:
        edges = np.linspace(float(np.min(combined)), float(np.max(combined)), 24)
    elif edges.size > 55:
        edges = np.linspace(float(np.min(combined)), float(np.max(combined)), 55)
    return edges


def make_pushforward_figure(
    freeze_file: str | Path = "freeze_samples.npz",
    *,
    output: str | Path = "pB_fig2_pushforward.pdf",
    write_current_alias: bool = True,
) -> list[Path]:
    sample = samples_from_npz(freeze_file)
    x_h = np.asarray(sample["XH"])
    h_phys = np.asarray(sample["H0P"])
    h_ladder = np.asarray(sample["H0L"])
    mu_x = float(sample["mu_x"])
    sigma_x = float(sample["sigma_x"])
    corr = float(sample["corr"])
    mu_p = float(sample["mu_phys"])
    sig_p = float(sample["sigma_phys"])
    mu_l = float(sample["mu_ladder"])
    sig_l = float(sample["sigma_ladder"])

    scale, exponent = _residual_scale(sigma_x)
    residual = (x_h - mu_x) * scale
    residual_sigma = sigma_x * scale

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(9.2, 4.15), layout="constrained"
    )

    # Actual posterior resamples.  The y-axis is only a linear change of units;
    # no jitter or synthetic width is introduced.
    ax_left.scatter(h_phys, residual, s=11, alpha=0.48, linewidths=0.0)
    ax_left.axhline(0.0, linewidth=0.8, linestyle="--")
    if residual_sigma > 0.0:
        ax_left.axhline(+residual_sigma, linewidth=0.65, linestyle=":")
        ax_left.axhline(-residual_sigma, linewidth=0.65, linestyle=":")
    ax_left.set_xlabel(r"$H_{0,\rm phys}\ [{\rm km\,s^{-1}\,Mpc^{-1}}]$")
    if exponent:
        ax_left.set_ylabel(
            rf"$10^{{{exponent}}}\,[X_H-\langle X_H\rangle]$"
        )
    else:
        ax_left.set_ylabel(r"$X_H-\langle X_H\rangle$")
    ax_left.grid(alpha=0.20)
    corr_text = "undefined" if not np.isfinite(corr) else f"{corr:+.3f}"
    ax_left.text(
        0.03,
        0.97,
        rf"actual posterior resamples" "\n"
        rf"$\langle X_H\rangle={mu_x:.7f}$" "\n"
        rf"$\sigma(X_H)={sigma_x:.2e}$" "\n"
        rf"${{\rm Corr}}(X_H,H_{{0,\rm phys}})={corr_text}$",
        transform=ax_left.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.82, edgecolor="0.7"),
    )

    # Empirical posterior densities from the same samples; no Gaussian model.
    edges = _common_histogram_edges(h_phys, h_ladder)
    ax_right.hist(
        h_phys,
        bins=edges,
        density=True,
        histtype="step",
        linewidth=1.8,
        label=rf"$H_{{0,\rm phys}}={mu_p:.2f}\pm{sig_p:.2f}$",
    )
    ax_right.hist(
        h_ladder,
        bins=edges,
        density=True,
        histtype="step",
        linewidth=1.8,
        label=rf"sample push-forward $={mu_l:.2f}\pm{sig_l:.2f}$",
    )
    # Rug marks show the underlying sample support without imposing a density
    # model.  Use a very short fixed height so the plot remains readable.
    rug_height = 0.018
    ax_right.vlines(h_phys, 0.0, rug_height, alpha=0.09, linewidth=0.45)
    ax_right.vlines(h_ladder, 0.0, rug_height, alpha=0.09, linewidth=0.45)
    ax_right.axvline(mu_p, linestyle="--", linewidth=0.8)
    ax_right.axvline(mu_l, linestyle="--", linewidth=0.8)
    ax_right.set_xlabel(r"$H_0\ [{\rm km\,s^{-1}\,Mpc^{-1}}]$")
    ax_right.set_ylabel("empirical posterior density")
    ax_right.set_ylim(bottom=0.0)
    ax_right.grid(alpha=0.20)
    ax_right.legend(frameon=False, fontsize=8.2, loc="upper center")
    p16, p50, p84 = _quantiles(h_ladder)
    ax_right.text(
        0.50,
        0.79,
        rf"push-forward $16/50/84\%={p16:.2f}/{p50:.2f}/{p84:.2f}$" "\n"
        r"No local-ladder data used",
        transform=ax_right.transAxes,
        ha="center",
        va="top",
        fontsize=8.2,
    )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    outputs = [output]
    fig.savefig(output, bbox_inches="tight")

    current_alias = output.with_name("pB_fig2_pushforward_v9.pdf")
    if write_current_alias and current_alias.resolve() != output.resolve():
        fig.savefig(current_alias, bbox_inches="tight")
        outputs.append(current_alias)
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("freeze_file", nargs="?", type=Path, default=Path("freeze_samples.npz"))
    parser.add_argument("--output", type=Path, default=Path("pB_fig2_pushforward.pdf"))
    parser.add_argument(
        "--no-v9-alias", action="store_true",
        help="do not also write pB_fig2_pushforward_v9.pdf",
    )
    args = parser.parse_args()
    outputs = make_pushforward_figure(
        args.freeze_file,
        output=args.output,
        write_current_alias=not args.no_v9_alias,
    )
    for out in outputs:
        print(f"saved {out}")


if __name__ == "__main__":
    main()
