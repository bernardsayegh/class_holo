#!/usr/bin/env python3
from pathlib import Path
import sys

pt = Path("source/perturbations.c")
t = pt.read_text()

if "SCR->CDM energy-transfer consistency (Option C)" in t:
    print("[SKIP] perturbations.c already patched")
    sys.exit(0)

# Find the marker line (your codebase uses this comment)
marker = "BEGIN HOLOGRAPHIC INTERACTION"
i = t.find(marker)
if i < 0:
    # fallback: search for any holographic block
    marker = "HOLOGRAPHIC INTERACTION"
    i = t.find(marker)
    if i < 0:
        print("[FATAL] Could not find HOLOGRAPHIC marker in source/perturbations.c")
        print("Run: grep -n \"HOLOGRAPHIC\" source/perturbations.c | head -50")
        sys.exit(1)

# Insert right AFTER the end of the comment line containing the marker
line_end = t.find("\n", i)
if line_end < 0:
    print("[FATAL] Unexpected file format")
    sys.exit(1)

insert = r'''
/* --- SCR->CDM energy-transfer consistency (Option C) -------------------------
   Uses pvecback[pba->index_bg_Q_scr_to_cdm_over_H] exported in background.c
   Frame: Q^μ || u_cdm (no momentum transfer in CDM frame).
   Default: δQ/Q = δ_cdm.
------------------------------------------------------------------------------- */
if (pba->has_super_schw_correction == _TRUE_ &&
    pba->index_bg_Q_scr_to_cdm_over_H >= 0 &&
    pba->super_schw_gamma != 0.) {

  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double rho_cdm = pvecback[pba->index_bg_rho_cdm];
  double Qscr_over_H = pvecback[pba->index_bg_Q_scr_to_cdm_over_H];

  if (rho_cdm > 0. && Qscr_over_H != 0.) {
    /* Qscr_over_H = Q/H  =>  aQ/rho = a*H*(Q/H)/rho */
    double coupling = a * H * Qscr_over_H / rho_cdm;

    double delta_cdm = y[pv->index_pt_delta_cdm];

    /* δQ/Q = δ_cdm (stable default) */
    double deltaQ_over_Q = delta_cdm;

    /* δ_c' += (aQ/rho_c) (δ_c - δQ/Q) */
    dy[pv->index_pt_delta_cdm] += coupling * (delta_cdm - deltaQ_over_Q);
  }
}
'''

t2 = t[:line_end+1] + insert + t[line_end+1:]
pt.write_text(t2)
print("[OK] perturbations.c patched")
