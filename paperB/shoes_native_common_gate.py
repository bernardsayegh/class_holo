#!/usr/bin/env python3
"""Native (M_B, a_B) common-gate test for the public SH0ES linear system.

The public SH0ES release represents the simultaneous distance-ladder fit as

    Y ~= theta @ L

with covariance C for Y.  This script solves the generalized least-squares
system, extracts M_B and 5 log10 H0, transforms to (M_B, a_B), and evaluates
both the common-gate amplitude and the normalization-independent null

    Delta M_B - 5 Delta a_B = 0.

It deliberately does not download data.  Supply the public Y/L/C products
from https://github.com/PantheonPlusSH0ES/DataRelease/SH0ES_Data .

Supported array formats: FITS primary image, NPY, NPZ, TXT, CSV.
No Astropy dependency is required for simple primary-image FITS files.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import struct
import sys
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy.linalg import cho_factor, cho_solve


class NativeLikelihoodError(RuntimeError):
    """Raised for malformed inputs or an ill-conditioned fit."""


def _parse_fits_value(raw: str) -> Any:
    raw = raw.split("/", 1)[0].strip()
    if not raw:
        return None
    if raw.startswith("'"):
        return raw.strip().strip("'").rstrip()
    if raw in {"T", "F"}:
        return raw == "T"
    try:
        return int(raw)
    except ValueError:
        try:
            return float(raw.replace("D", "E"))
        except ValueError:
            return raw


def read_primary_fits(path: pathlib.Path) -> np.ndarray:
    """Read a simple FITS primary image without external dependencies."""
    with path.open("rb") as f:
        cards: list[str] = []
        found_end = False
        while not found_end:
            block = f.read(2880)
            if len(block) != 2880:
                raise NativeLikelihoodError(f"Truncated FITS header: {path}")
            for i in range(0, 2880, 80):
                card = block[i : i + 80].decode("ascii", errors="replace")
                cards.append(card)
                if card.startswith("END"):
                    found_end = True
                    break

        header: dict[str, Any] = {}
        for card in cards:
            key = card[:8].strip()
            if key == "END":
                break
            if card[8:10] == "= ":
                header[key] = _parse_fits_value(card[10:])

        bitpix = int(header.get("BITPIX", 0))
        naxis = int(header.get("NAXIS", 0))
        if naxis < 1:
            raise NativeLikelihoodError(f"FITS primary HDU has no image: {path}")
        axes = [int(header[f"NAXIS{i}"]) for i in range(1, naxis + 1)]
        count = int(np.prod(axes, dtype=np.int64))
        dtype_map = {
            8: np.dtype(">u1"),
            16: np.dtype(">i2"),
            32: np.dtype(">i4"),
            64: np.dtype(">i8"),
            -32: np.dtype(">f4"),
            -64: np.dtype(">f8"),
        }
        if bitpix not in dtype_map:
            raise NativeLikelihoodError(f"Unsupported BITPIX={bitpix} in {path}")
        data = np.fromfile(f, dtype=dtype_map[bitpix], count=count)
        if data.size != count:
            raise NativeLikelihoodError(f"Truncated FITS data array: {path}")
        # FITS NAXIS1 is the fastest-varying axis.
        arr = data.reshape(tuple(reversed(axes))).astype(float, copy=False)
        bscale = float(header.get("BSCALE", 1.0))
        bzero = float(header.get("BZERO", 0.0))
        return arr * bscale + bzero


def load_array(path_like: str) -> np.ndarray:
    path = pathlib.Path(path_like)
    if not path.exists():
        raise NativeLikelihoodError(f"File not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".fits", ".fit", ".fts"}:
        return read_primary_fits(path)
    if suffix == ".npy":
        return np.load(path)
    if suffix == ".npz":
        archive = np.load(path)
        if len(archive.files) != 1:
            raise NativeLikelihoodError(
                f"NPZ must contain exactly one array; found {archive.files}"
            )
        return archive[archive.files[0]]
    delimiter = "," if suffix == ".csv" else None
    return np.loadtxt(path, delimiter=delimiter)


def _as_vector(a: np.ndarray, name: str) -> np.ndarray:
    v = np.asarray(a, dtype=float).squeeze()
    if v.ndim != 1:
        raise NativeLikelihoodError(f"{name} must be a vector; got shape {a.shape}")
    return v


def _as_cov(a: np.ndarray, n: int, name: str) -> np.ndarray:
    c = np.asarray(a, dtype=float).squeeze()
    if c.ndim == 1 and c.size == n * n:
        c = c.reshape(n, n)
    if c.shape != (n, n):
        raise NativeLikelihoodError(f"{name} must have shape {(n, n)}; got {c.shape}")
    return 0.5 * (c + c.T)


@dataclass
class GLSResult:
    theta: np.ndarray
    covariance: np.ndarray
    chi2: float
    dof: int
    condition_fisher: float


def solve_shoes_gls(Y: np.ndarray, L: np.ndarray, C: np.ndarray) -> GLSResult:
    """Solve Y = theta @ L + noise using stable Cholesky GLS."""
    y = _as_vector(Y, "Y")
    design = np.asarray(L, dtype=float).squeeze()
    if design.ndim != 2:
        raise NativeLikelihoodError(f"L must be a matrix; got shape {design.shape}")
    # Public SH0ES convention is L=(n_parameter,n_datum), theta @ L.
    if design.shape[1] != y.size and design.shape[0] == y.size:
        design = design.T
    if design.shape[1] != y.size:
        raise NativeLikelihoodError(
            f"L is incompatible with Y: L={design.shape}, Y={y.shape}"
        )
    cov_y = _as_cov(C, y.size, "C")
    try:
        cf = cho_factor(cov_y, lower=True, check_finite=True)
    except Exception as exc:
        raise NativeLikelihoodError(f"C is not positive definite: {exc}") from exc

    cinv_y = cho_solve(cf, y)
    cinv_lt = cho_solve(cf, design.T)
    fisher = design @ cinv_lt
    rhs = design @ cinv_y
    cond = float(np.linalg.cond(fisher))
    if not np.isfinite(cond) or cond > 1e16:
        raise NativeLikelihoodError(f"Ill-conditioned Fisher matrix: cond={cond:.3e}")
    theta = np.linalg.solve(fisher, rhs)
    cov_theta = np.linalg.inv(fisher)
    residual = y - theta @ design
    chi2 = float(residual @ cho_solve(cf, residual))
    dof = int(y.size - theta.size)
    return GLSResult(theta, cov_theta, chi2, dof, cond)


def native_coordinates(
    theta: np.ndarray, cov_theta: np.ndarray, i_mb: int, i_h5: int
) -> tuple[np.ndarray, np.ndarray]:
    """Transform (M_B, 5 log10 H0) to (M_B, a_B)."""
    n = theta.size
    for index, name in [(i_mb, "M_B"), (i_h5, "5log10H0")]:
        if not 0 <= index < n:
            raise NativeLikelihoodError(f"Index {index} for {name} outside 0..{n-1}")
    mb = float(theta[i_mb])
    h5 = float(theta[i_h5])
    ab = (h5 - mb - 25.0) / 5.0
    mean = np.array([mb, ab], dtype=float)
    sub = cov_theta[np.ix_([i_mb, i_h5], [i_mb, i_h5])]
    jac = np.array([[1.0, 0.0], [-0.2, 0.2]])
    cov = jac @ sub @ jac.T
    return mean, cov


def read_physical_reference(path: str) -> tuple[np.ndarray, np.ndarray]:
    payload = json.loads(pathlib.Path(path).read_text())
    mean = np.asarray(payload["mean"], dtype=float)
    cov = np.asarray(payload["covariance"], dtype=float)
    if mean.shape != (2,) or cov.shape != (2, 2):
        raise NativeLikelihoodError(
            "Physical reference JSON needs mean=[M_B,a_B] and 2x2 covariance"
        )
    return mean, 0.5 * (cov + cov.T)


@dataclass
class CommonGateResult:
    delta: list[float]
    covariance_delta: list[list[float]]
    x_hat: float
    sigma_x: float
    a_map_hat: float | None
    delta_mb_minus_5_delta_ab: float
    sigma_null: float
    chi2_null: float
    pull_null: float
    chi2_common_mode: float
    predicted_delta_for_x: list[float] | None


def fit_common_gate(
    observed_mean: np.ndarray,
    observed_cov: np.ndarray,
    physical_mean: np.ndarray,
    physical_cov: np.ndarray,
    x_prediction: float | None,
) -> CommonGateResult:
    delta = observed_mean - physical_mean
    cov_delta = observed_cov + physical_cov
    cf = cho_factor(cov_delta, lower=True, check_finite=True)
    # z=(M_B,a_B), so z shift per unit X is (5/ln10,1/ln10).
    s = np.array([5.0 / math.log(10.0), 1.0 / math.log(10.0)])
    cinv_s = cho_solve(cf, s)
    denom = float(s @ cinv_s)
    x_hat = float(s @ cho_solve(cf, delta) / denom)
    sigma_x = math.sqrt(1.0 / denom)
    resid = delta - s * x_hat
    chi2_common = float(resid @ cho_solve(cf, resid))

    nvec = np.array([1.0, -5.0])
    dnull = float(nvec @ delta)
    varnull = float(nvec @ cov_delta @ nvec)
    if varnull <= 0:
        raise NativeLikelihoodError("Non-positive null-mode variance")
    sigma_null = math.sqrt(varnull)
    pull = dnull / sigma_null
    pred = None if x_prediction is None else (s * x_prediction).tolist()
    amap = None if x_prediction in {None, 0.0} else 2.0 * x_hat / x_prediction
    return CommonGateResult(
        delta=delta.tolist(),
        covariance_delta=cov_delta.tolist(),
        x_hat=x_hat,
        sigma_x=sigma_x,
        a_map_hat=amap,
        delta_mb_minus_5_delta_ab=dnull,
        sigma_null=sigma_null,
        chi2_null=pull * pull,
        pull_null=pull,
        chi2_common_mode=chi2_common,
        predicted_delta_for_x=pred,
    )


def self_test() -> dict[str, Any]:
    rng = np.random.default_rng(240713)
    npar, ndata = 7, 80
    L = rng.normal(size=(npar, ndata))
    # Ensure excellent rank.
    theta_true = rng.normal(size=npar)
    A = rng.normal(size=(ndata, ndata))
    C = (A @ A.T) / ndata + np.eye(ndata) * 0.2
    noise = rng.multivariate_normal(np.zeros(ndata), C)
    Y = theta_true @ L + noise
    result = solve_shoes_gls(Y, L, C)
    # Verify the normal equations, not exact truth recovery in one noisy sample.
    normal = L @ np.linalg.solve(C, Y - result.theta @ L)
    normal_error = float(np.max(np.abs(normal)))
    if normal_error > 1e-8:
        raise NativeLikelihoodError(f"Self-test normal equations failed: {normal_error}")

    x_true = 0.0359
    s = np.array([5 / math.log(10), 1 / math.log(10)])
    phys_mean = np.array([-19.3, 0.71])
    obs_mean = phys_mean + s * x_true
    obs_cov = np.diag([0.02**2, 0.004**2])
    phys_cov = np.diag([0.01**2, 0.002**2])
    cg = fit_common_gate(obs_mean, obs_cov, phys_mean, phys_cov, x_true)
    if abs(cg.x_hat - x_true) > 1e-12 or abs(cg.delta_mb_minus_5_delta_ab) > 1e-12:
        raise NativeLikelihoodError("Self-test common-gate recovery failed")
    return {
        "status": "passed",
        "normal_equation_max_abs": normal_error,
        "synthetic_x_true": x_true,
        "synthetic_x_hat": cg.x_hat,
        "synthetic_null": cg.delta_mb_minus_5_delta_ab,
    }


def cmd_fit(args: argparse.Namespace) -> int:
    Y = load_array(args.y)
    L = load_array(args.l)
    C = load_array(args.c)
    gls = solve_shoes_gls(Y, L, C)
    native_mean, native_cov = native_coordinates(
        gls.theta, gls.covariance, args.i_mb, args.i_h5
    )
    output: dict[str, Any] = {
        "gls": {
            "chi2": gls.chi2,
            "dof": gls.dof,
            "condition_fisher": gls.condition_fisher,
            "n_parameter": int(gls.theta.size),
            "native_parameter_indices": {"M_B": args.i_mb, "5log10H0": args.i_h5},
        },
        "native_observed": {
            "mean_M_B_a_B": native_mean.tolist(),
            "covariance_M_B_a_B": native_cov.tolist(),
            "sigma_M_B_a_B": np.sqrt(np.diag(native_cov)).tolist(),
        },
    }
    if args.physical_reference:
        phys_mean, phys_cov = read_physical_reference(args.physical_reference)
        cg = fit_common_gate(
            native_mean, native_cov, phys_mean, phys_cov, args.x_prediction
        )
        output["physical_reference"] = {
            "mean_M_B_a_B": phys_mean.tolist(),
            "covariance_M_B_a_B": phys_cov.tolist(),
        }
        output["common_gate"] = asdict(cg)
    else:
        output["note"] = (
            "No physical reference supplied. The public ladder coordinates were fitted, "
            "but common-gate shifts require an independently frozen non-ladder reference."
        )
    text = json.dumps(output, indent=2, sort_keys=True)
    if args.output:
        pathlib.Path(args.output).write_text(text + "\n")
    print(text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fit = sub.add_parser("fit", help="Fit the public SH0ES linear system")
    fit.add_argument("--y", required=True, help="SH0ES Y vector")
    fit.add_argument("--l", required=True, help="SH0ES L/design matrix")
    fit.add_argument("--c", required=True, help="SH0ES covariance matrix")
    fit.add_argument("--i-mb", type=int, default=42, help="zero-based index of M_B")
    fit.add_argument(
        "--i-h5", type=int, default=46, help="zero-based index of 5 log10 H0"
    )
    fit.add_argument(
        "--physical-reference",
        help="JSON containing mean=[M_B,a_B] and covariance=[[...],[...]]",
    )
    fit.add_argument(
        "--x-prediction", type=float, help="frozen non-ladder X_H prediction"
    )
    fit.add_argument("--output", help="write JSON result to this path")
    fit.set_defaults(func=cmd_fit)

    test = sub.add_parser("self-test", help="Run synthetic algebraic tests")
    test.set_defaults(func=lambda _: (print(json.dumps(self_test(), indent=2)) or 0))
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except NativeLikelihoodError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
