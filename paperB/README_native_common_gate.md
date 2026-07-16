# Native SH0ES common-gate likelihood module

`shoes_native_common_gate.py` implements the native two-coordinate test used in Paper B.
It solves the public SH0ES generalized least-squares system

\[
Y\simeq \theta L,
\]

extracts \((M_B,5\log_{10}H_0)\), transforms to

\[
a_B=\frac{5\log_{10}H_0-M_B-25}{5},
\]

and evaluates the common-gate shift

\[
\delta z = X_H\left(\frac{5}{\ln 10},\frac{1}{\ln 10}\right),
\qquad z=(M_B,a_B).
\]

The normalization-independent null is

\[
D_{\rm com}\equiv\delta M_B-5\delta a_B=0.
\]

## Inputs

Use the public data products in the official Pantheon+SH0ES repository:

- `Y`: data vector;
- `L`: design matrix in the public convention `theta @ L`;
- `C`: covariance matrix.

The script accepts FITS primary images, NPY, NPZ, TXT, and CSV. A minimal FITS reader is bundled, so Astropy is not required for the simple primary-image products.

The default zero-based parameter indices are

- `M_B`: 42;
- `5 log10 H0`: 46.

Verify these against the release version used in an analysis. They can be changed with `--i-mb` and `--i-h5`.

## Physical reference

The common-gate amplitude requires an independently frozen non-ladder reference in native coordinates. Supply JSON of the form

```json
{
  "mean": [-19.30, 0.71],
  "covariance": [[0.0004, 0.0], [0.0, 0.000016]]
}
```

The example numbers above are placeholders only. They are not Paper-B results.

## Commands

Synthetic algebraic test:

```bash
python shoes_native_common_gate.py self-test
```

Fit the public system and report native coordinates:

```bash
python shoes_native_common_gate.py fit \
  --y SH0ES_Y.fits \
  --l SH0ES_L.fits \
  --c SH0ES_C.fits
```

Fit and test the frozen common-gate prediction:

```bash
python shoes_native_common_gate.py fit \
  --y SH0ES_Y.fits \
  --l SH0ES_L.fits \
  --c SH0ES_C.fits \
  --physical-reference nonladder_native_reference.json \
  --x-prediction 0.0359 \
  --output native_common_gate_result.json
```

## Outputs

The JSON output includes:

- generalized least-squares \(\chi^2\), degrees of freedom, and Fisher condition number;
- fitted \((M_B,a_B)\) mean and covariance;
- fitted common-gate amplitude \(\widehat X\) and uncertainty;
- the native null \(D_{\rm com}\), its uncertainty, pull, and \(\chi^2\);
- the residual \(\chi^2\) orthogonal to the common mode;
- optionally, the implied mapping amplitude relative to a frozen \(X_H\).

## Scope

This package completes the statistical implementation and passes a synthetic recovery test. It does not ship or fabricate the public SH0ES covariance products, and it does not substitute placeholder values for a real analysis. A publication result requires the official data release and a separately frozen non-ladder native reference.
