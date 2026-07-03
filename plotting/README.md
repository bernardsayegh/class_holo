# plotting/ -- paper figure scripts (v5 grid)

Script names match the \includegraphics filenames in the paper: each
script X.py produces X.pdf (+ .png) and prints a summary block for
caption verification. Chain-dependent scripts resolve roots from
~/class_holo_test/chains_test automatically (v5 naming, glob fallback).

| script | paper figure | needs |
|---|---|---|
| posteriors.py | fig:posteriors | LCDM, A, B(A2) chains |
| S8_barchart.py | fig:S8_barchart | LCDM, A, B(A2) chains |
| contours_beta_S8.py | fig:S8_contours | A(beta-free) chain |
| triangle_S8_Om_H0.py | fig:triangle | LCDM, A, B(A2) chains |
| S8_fclust.py | fig:S8_panels | CLASS build + A(fc-free) chain |
| spectra_Pk.py | fig:spectra_Pk | CLASS build (fixed cosmology) |
| spectra_CMB.py | fig:spectra_CMB | CLASS build (theta_s-matched) |
| fsigma8_comparison.py | fig:fsigma8 | CLASS build |
| evolution_density.py | fig:evolution | CLASS build |
| desi_bao_dr2.py | fig:desi_bao | CLASS build |
| contours_H0_S8.py | (supplementary) | LCDM, A, B(A2) chains |
| holographic_data_bracket_expanded.py | (supplementary) | none |

Not yet included: sweep_rate_omega.py and evolution_Sz.py (analytic
figures fig:sweep_rate and fig:Sz). CLASS-based scripts must run from
the repo root (they add ./python to sys.path). Reference bands: SH0ES
R25 (73.49 +/- 0.93), DES-Y3, growth fsigma8 S8 = 0.747 +/- 0.029
(Nguyen+ 2023), Planck and Pantheon+ Omega_m. LCDM line conventions:
0.834 = Planck-only, 0.787 = anchored joint fit.
