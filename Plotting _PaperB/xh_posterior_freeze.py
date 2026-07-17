#!/usr/bin/env python3
"""Freeze the non-ladder posterior and reconstruct the Paper-B exposure.

The script keeps its original filename and command-line entry point, but fixes
several issues in the earlier version:

* posterior sample weights are respected through deterministic systematic
  resampling;
* the sweep kernel uses the intended two-fluid fractions
  ``rho_lambda/(rho_lambda+rho_b+rho_cdm)``;
* the production smooth gate is reproduced by default;
* redshift interpolation is always performed on an increasing grid;
* the near-zero-width ``X_H`` scatter is no longer used for Figure 2;
* a JSON summary and SHA-256 chain manifest are written for provenance.

Example
-------

    python3 xh_posterior_freeze.py chains_freeze/modelB_noshoes.1.txt

Outputs (in ``--output-dir``):

* ``freeze_samples.npz``
* ``freeze_summary.json``
* ``freeze_chain_manifest.sha256``
* ``pB_fig1_capacity_history.pdf``
* ``pB_fig2_pushforward.pdf`` and ``pB_fig2_pushforward_v9.pdf``
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_LIGHT_KM_S = 299792.458


@dataclass(frozen=True)
class FreezeConfig:
    burn_in: float = 0.30
    n_samples: int = 401
    seed: int = 1729
    beta: float = 1.0 / 12.0
    delta_s: float = 0.01
    hard_gate: bool = False
    z_plot_max: float = 3.0


@dataclass(frozen=True)
class ExposureCurve:
    z: np.ndarray          # increasing redshift
    sweep: np.ndarray
    excess: np.ndarray
    cumulative: np.ndarray # X accumulated from early times to epoch z


@dataclass(frozen=True)
class ExposureRecord:
    x_h: float
    z_entry: float
    z_exit: float
    h0_phys: float
    curve: ExposureCurve


def normalize_chain_root(path: str | Path) -> str:
    text = str(path)
    match = re.match(r"^(.*)\.([0-9]+)\.txt$", text)
    if match:
        return match.group(1)
    if text.endswith(".txt"):
        return text[:-4]
    return text


def discover_chain_files(root_or_file: str | Path) -> list[Path]:
    root = normalize_chain_root(root_or_file)
    files = [Path(p) for p in sorted(glob.glob(root + ".[0-9]*.txt"))]
    if not files:
        direct = Path(root_or_file)
        if direct.is_file():
            files = [direct]
    if not files:
        raise FileNotFoundError(f"no chain files match {root}.N.txt")
    return files


def read_chain_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline()
    if not first.lstrip().startswith("#"):
        raise ValueError(f"chain file has no commented header: {path}")
    names = first.lstrip("#").split()
    if not names:
        raise ValueError(f"empty chain header: {path}")
    return names


def load_post_burnin(files: Iterable[Path], burn_in: float) -> tuple[list[str], np.ndarray]:
    if not 0.0 <= burn_in < 1.0:
        raise ValueError("burn_in must lie in [0,1)")
    files = list(files)
    names = read_chain_header(files[0])
    parts: list[np.ndarray] = []
    for path in files:
        if read_chain_header(path) != names:
            raise ValueError(f"chain header differs in {path}")
        data = np.loadtxt(path, comments="#", ndmin=2)
        if data.shape[1] != len(names):
            raise ValueError(
                f"{path}: {data.shape[1]} columns but header has {len(names)} names"
            )
        start = int(np.floor(burn_in * data.shape[0]))
        kept = data[start:]
        if kept.size:
            parts.append(kept)
    if not parts:
        raise ValueError("no samples remain after burn-in")
    return names, np.vstack(parts)


def _find_column(columns: dict[str, int], aliases: Iterable[str], *, required: bool = True) -> int | None:
    for name in aliases:
        if name in columns:
            return columns[name]
    if required:
        raise KeyError(f"none of the required columns are present: {list(aliases)}")
    return None


def posterior_weights(raw: np.ndarray, columns: dict[str, int]) -> np.ndarray:
    index = _find_column(columns, ("weight", "weights", "sample_weight"), required=False)
    if index is None:
        weights = np.ones(raw.shape[0], dtype=float)
    else:
        weights = np.asarray(raw[:, index], dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError("posterior weights have zero total")
    return weights / total


def systematic_resample(probabilities: np.ndarray, n_samples: int, seed: int) -> np.ndarray:
    if n_samples < 2:
        raise ValueError("n_samples must be at least two")
    probabilities = np.asarray(probabilities, dtype=float)
    probabilities = probabilities / probabilities.sum()
    rng = np.random.default_rng(seed)
    start = rng.random() / n_samples
    positions = start + np.arange(n_samples) / n_samples
    cumulative = np.cumsum(probabilities)
    cumulative[-1] = 1.0
    return np.searchsorted(cumulative, positions, side="right")


def _parameter(row: np.ndarray, columns: dict[str, int], aliases: Iterable[str]) -> float:
    index = _find_column(columns, aliases)
    value = float(row[index])
    if not np.isfinite(value):
        raise ValueError(f"non-finite chain parameter in {list(aliases)}")
    return value


def row_to_class_parameters(row: np.ndarray, columns: dict[str, int]) -> dict[str, float]:
    h0_index = _find_column(columns, ("H0", "H_0"), required=False)
    if h0_index is not None:
        h0 = float(row[h0_index])
    else:
        h0 = 100.0 * _parameter(row, columns, ("h",))

    a_s_index = _find_column(columns, ("A_s", "As"), required=False)
    if a_s_index is not None:
        a_s = float(row[a_s_index])
    else:
        log_a = _parameter(row, columns, ("logA", "ln10^{10}A_s", "ln10^10A_s"))
        a_s = 1.0e-10 * np.exp(log_a)

    return {
        "omega_b": _parameter(row, columns, ("omega_b", "ombh2")),
        "omega_cdm": _parameter(row, columns, ("omega_cdm", "omch2")),
        "H0": h0,
        "tau_reio": _parameter(row, columns, ("tau_reio", "tau")),
        "n_s": _parameter(row, columns, ("n_s", "ns")),
        "A_s": a_s,
    }


def _background_column(background: dict[str, np.ndarray], *names: str) -> np.ndarray:
    for name in names:
        if name in background:
            return np.asarray(background[name], dtype=float)
    raise KeyError(f"none of the CLASS background columns are present: {names}")


def _regulated_excess(sweep: np.ndarray, delta_s: float, hard_gate: bool) -> np.ndarray:
    safe = np.maximum(sweep, 1.0e-15)
    hard = np.where(sweep > 1.0, 1.0 - 1.0 / safe, 0.0)
    if hard_gate:
        return hard
    if delta_s <= 0.0:
        raise ValueError("delta_s must be positive for the smooth gate")
    gate = 0.5 * (1.0 + np.tanh((sweep - 1.0) / delta_s))
    return gate * hard


def _crossings(z: np.ndarray, values: np.ndarray, level: float = 1.0) -> list[float]:
    z = np.asarray(z, dtype=float)
    y = np.asarray(values, dtype=float) - level
    order = np.argsort(z)
    z, y = z[order], y[order]
    out: list[float] = []
    for i in np.where(y[:-1] * y[1:] <= 0.0)[0]:
        if y[i] == y[i + 1]:
            continue
        frac = -y[i] / (y[i + 1] - y[i])
        if 0.0 <= frac <= 1.0:
            out.append(float(z[i] + frac * (z[i + 1] - z[i])))
    # Remove near-duplicates caused by exact grid hits.
    unique: list[float] = []
    for value in sorted(out):
        if not unique or abs(value - unique[-1]) > 1.0e-5:
            unique.append(value)
    return unique


def evaluate_sample(parameters: dict[str, float], config: FreezeConfig) -> ExposureRecord:
    try:
        from classy import Class
    except ImportError as exc:
        raise RuntimeError(
            "classy is required for chain reconstruction; run inside the modified CLASS environment"
        ) from exc

    cosmo = Class()
    settings: dict[str, float | int | str] = {
        **parameters,
        "YHe": 0.2454,
        "N_ur": 2.0328,
        "N_ncdm": 1,
        "m_ncdm": 0.06,
        "interaction_beta": config.beta,
        "interaction_ieff_type": 4,
        "f_clust": 0.0,
        "super_schwarzschild_correction": "yes",
        "super_schw_amp": 0.0,
        "super_schw_Amap": 2.0,
        "super_schw_deltaS": config.delta_s,
        "super_schw_gamma": 2.0,
        "super_schw_ode": 0,
        "output": "",
    }
    try:
        cosmo.set(settings)
        cosmo.compute()
        bg = cosmo.get_background()
        z_raw = _background_column(bg, "z")
        rho_l = _background_column(bg, "(.)rho_lambda", "(.)rho_fld")
        rho_b = _background_column(bg, "(.)rho_b")
        rho_c = _background_column(bg, "(.)rho_cdm")
        denom = rho_l + rho_b + rho_c
        omega_l_2f = np.divide(rho_l, denom, out=np.zeros_like(rho_l), where=denom > 0.0)
        sweep_raw = 4.5 * omega_l_2f * (1.0 - omega_l_2f)
        excess_raw = _regulated_excess(sweep_raw, config.delta_s, config.hard_gate)

        a_raw = 1.0 / (1.0 + z_raw)
        lna_raw = np.log(a_raw)
        time_order = np.argsort(lna_raw)  # early -> late
        lna_t = lna_raw[time_order]
        sweep_t = sweep_raw[time_order]
        excess_t = excess_raw[time_order]
        z_t = z_raw[time_order]
        cumulative_t = np.concatenate(
            [[0.0], np.cumsum(0.5 * (excess_t[1:] + excess_t[:-1]) * np.diff(lna_t))]
        )
        x_h = float(cumulative_t[-1])

        # Interpolation routines require increasing redshift.  Reversing the
        # early->late arrays gives z=0 first and keeps cumulative(z=0)=X_H.
        z_order = np.argsort(z_t)
        curve = ExposureCurve(
            z=np.asarray(z_t[z_order], dtype=float),
            sweep=np.asarray(sweep_t[z_order], dtype=float),
            excess=np.asarray(excess_t[z_order], dtype=float),
            cumulative=np.asarray(cumulative_t[z_order], dtype=float),
        )
        roots = _crossings(curve.z, curve.sweep, 1.0)
        positive = [value for value in roots if value >= 0.0]
        z_exit = min(positive) if positive else float("nan")
        z_entry = max(positive) if positive else float("nan")
        h0_phys = float(cosmo.Hubble(0.0) * C_LIGHT_KM_S)
        return ExposureRecord(x_h=x_h, z_entry=z_entry, z_exit=z_exit,
                              h0_phys=h0_phys, curve=curve)
    finally:
        try:
            cosmo.struct_cleanup()
            cosmo.empty()
        except Exception:
            pass


def evaluate_chain(
    root_or_file: str | Path,
    config: FreezeConfig,
) -> tuple[list[Path], list[str], np.ndarray, np.ndarray, list[ExposureRecord]]:
    files = discover_chain_files(root_or_file)
    names, raw = load_post_burnin(files, config.burn_in)
    columns = {name: i for i, name in enumerate(names)}
    probabilities = posterior_weights(raw, columns)
    indices = systematic_resample(probabilities, config.n_samples, config.seed)
    rows = raw[indices]

    records: list[ExposureRecord] = []
    for number, row in enumerate(rows, start=1):
        records.append(evaluate_sample(row_to_class_parameters(row, columns), config))
        if number == 1 or number % 25 == 0 or number == len(rows):
            print(f"reconstructed {number}/{len(rows)} posterior samples", flush=True)
    return files, names, rows, raw, records


def _summary(records: list[ExposureRecord]) -> dict[str, object]:
    x_h = np.array([r.x_h for r in records])
    z_entry = np.array([r.z_entry for r in records])
    z_exit = np.array([r.z_exit for r in records])
    h0_phys = np.array([r.h0_phys for r in records])
    h0_ladder = h0_phys * np.exp(2.0 * x_h)

    def stats(array: np.ndarray) -> dict[str, float]:
        finite = array[np.isfinite(array)]
        return {
            "mean": float(np.mean(finite)),
            "std": float(np.std(finite, ddof=1)),
            "p16": float(np.percentile(finite, 16.0)),
            "median": float(np.percentile(finite, 50.0)),
            "p84": float(np.percentile(finite, 84.0)),
        }

    corr = float(np.corrcoef(x_h, h0_phys)[0, 1]) if np.std(x_h) > 0.0 else float("nan")
    return {
        "n_samples": len(records),
        "XH": stats(x_h),
        "z_entry": stats(z_entry),
        "z_exit": stats(z_exit),
        "H0_phys": stats(h0_phys),
        "H0_ladder": stats(h0_ladder),
        "corr_XH_H0_phys": corr,
        "mean_mapping_factor": float(np.mean(np.exp(2.0 * x_h))),
    }


def write_manifest(files: Iterable[Path], output: Path) -> None:
    lines: list[str] = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.resolve()}")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_capacity_history(
    records: list[ExposureRecord],
    *,
    output: str | Path = "pB_fig1_capacity_history.pdf",
    z_max: float = 3.0,
) -> Path:
    if not records:
        raise ValueError("at least one exposure record is required")
    z_grid = np.linspace(0.0, z_max, 500)

    def matrix(attribute: str) -> np.ndarray:
        values = []
        for record in records:
            curve = record.curve
            values.append(np.interp(z_grid, curve.z, getattr(curve, attribute)))
        return np.asarray(values)

    sweep = matrix("sweep")
    excess = matrix("excess")
    cumulative = matrix("cumulative")
    z_entry = np.array([r.z_entry for r in records])
    z_exit = np.array([r.z_exit for r in records])

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6.0, 7.0), layout="constrained")
    for ax, values, label in zip(
        axes,
        (sweep, excess, cumulative),
        (r"$\mathcal{S}(z)$", r"$\Delta(z)$", r"$X_H^{\rm cum}(z)$"),
    ):
        lo, median, hi = np.percentile(values, [16.0, 50.0, 84.0], axis=0)
        ax.fill_between(z_grid, lo, hi, alpha=0.28)
        ax.plot(z_grid, median, linewidth=1.6)
        ax.set_ylabel(label)
        ax.grid(alpha=0.18)
    axes[0].axhline(1.0, linestyle=":", linewidth=0.8)
    for ax in axes:
        if np.any(np.isfinite(z_entry)):
            ax.axvline(np.nanmedian(z_entry), linestyle="--", linewidth=0.75)
        if np.any(np.isfinite(z_exit)):
            ax.axvline(np.nanmedian(z_exit), linestyle="--", linewidth=0.75)
    axes[2].set_xlabel(r"redshift $z$")
    axes[2].set_ylim(bottom=0.0)
    axes[2].invert_xaxis()

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def save_freeze(
    records: list[ExposureRecord],
    rows: np.ndarray,
    names: list[str],
    *,
    output_dir: Path,
    config: FreezeConfig,
) -> tuple[Path, Path]:
    x_h = np.array([r.x_h for r in records])
    z_entry = np.array([r.z_entry for r in records])
    z_exit = np.array([r.z_exit for r in records])
    h0_phys = np.array([r.h0_phys for r in records])
    npz_path = output_dir / "freeze_samples.npz"
    np.savez(
        npz_path,
        XH=x_h,
        H0P=h0_phys,
        ZIN=z_entry,
        ZOUT=z_exit,
        chain_rows=rows,
        chain_columns=np.asarray(names, dtype="U"),
        beta=np.array(config.beta),
        delta_s=np.array(config.delta_s),
        hard_gate=np.array(config.hard_gate),
        burn_in=np.array(config.burn_in),
        seed=np.array(config.seed),
    )
    summary = _summary(records)
    summary.update({
        "configuration": {
            "burn_in": config.burn_in,
            "n_samples": config.n_samples,
            "seed": config.seed,
            "beta": config.beta,
            "delta_s": config.delta_s,
            "hard_gate": config.hard_gate,
        }
    })
    json_path = output_dir / "freeze_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    return npz_path, json_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain", type=Path, help="chain root or any root.N.txt file")
    parser.add_argument("--burn-in", type=float, default=0.30)
    parser.add_argument("--samples", type=int, default=401)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--beta", type=float, default=1.0 / 12.0)
    parser.add_argument("--delta-s", type=float, default=0.01)
    parser.add_argument("--hard-gate", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    config = FreezeConfig(
        burn_in=args.burn_in,
        n_samples=args.samples,
        seed=args.seed,
        beta=args.beta,
        delta_s=args.delta_s,
        hard_gate=args.hard_gate,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files, names, rows, raw, records = evaluate_chain(args.chain, config)
    npz_path, summary_path = save_freeze(records, rows, names, output_dir=args.output_dir, config=config)
    manifest_path = args.output_dir / "freeze_chain_manifest.sha256"
    write_manifest(files, manifest_path)

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(json.dumps(summary, indent=2))

    columns = {name: i for i, name in enumerate(names)}
    stored_index = _find_column(columns, ("H0_local", "H0_ladder"), required=False)
    if stored_index is not None:
        stored = rows[:, stored_index]
        predicted = np.array([r.h0_phys * np.exp(2.0 * r.x_h) for r in records])
        finite = np.isfinite(stored) & (stored != 0.0)
        relative = np.abs(predicted[finite] - stored[finite]) / np.abs(stored[finite])
        print(
            "cross-check against stored mapped H0: "
            f"max={np.max(relative):.3e}, mean={np.mean(relative):.3e}"
        )

    if not args.no_figures:
        figure1 = plot_capacity_history(
            records,
            output=args.output_dir / "pB_fig1_capacity_history.pdf",
            z_max=config.z_plot_max,
        )
        from pB_fig2_pushforward import make_pushforward_figure
        figure2 = make_pushforward_figure(
            npz_path,
            output=args.output_dir / "pB_fig2_pushforward.pdf",
            write_current_alias=True,
        )
        print(f"saved {figure1}")
        for path in figure2:
            print(f"saved {path}")
    print(f"saved {npz_path}")
    print(f"saved {summary_path}")
    print(f"saved {manifest_path}")


if __name__ == "__main__":
    main()
