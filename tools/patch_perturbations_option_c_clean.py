#!/usr/bin/env python3
from pathlib import Path
import sys

p = Path("source/perturbations.c")
lines = p.read_text().splitlines(True)

# Find the holographic perturbations header
hdr = None
for i, ln in enumerate(lines):
    if "HOLOGRAPHIC INTERACTION (PERTURBATIONS)" in ln:
        hdr = i
        break
if hdr is None:
    print("[FATAL] Can't find 'HOLOGRAPHIC INTERACTION (PERTURBATIONS)'")
    sys.exit(1)

# Identify the "doc-star" block: consecutive lines after header that start with whitespace + '*'
doc_start = None
for j in range(hdr + 1, len(lines)):
    s = lines[j].lstrip()
    if s.startswith("*"):
        doc_start = j
        break
    # if we hit real code before any star lines, no doc block
    if s and not s.startswith(("/*", "//")):
        break

doc_end = None
insert_after = hdr + 1

if doc_start is not None:
    # Walk until star-lines stop (allow blank lines)
    k = doc_start
    while k < len(lines):
        s = lines[k].lstrip()
        if s.startswith("*") or s == "" or s == "\n":
            k += 1
            continue
        break
    doc_end = k  # first non-doc line index

    # Ensure doc block is wrapped in /* ... */
    # Insert opening /* right before first '*' line if not already inside a block comment
    # (If previous line already has /* without closing, we won't add)
    prev_chunk = "".join(lines[max(hdr, doc_start-5):doc_start])
    if "/*" not in prev_chunk or prev_chunk.rfind("/*") < prev_chunk.rfind("*/"):
        lines.insert(doc_start, "/*\n")
        doc_end += 1

    # Insert closing */ after doc block if not already closed inside
    chunk = "".join(lines[doc_start:doc_end])
    if "*/" not in chunk:
        lines.insert(doc_end, "*/\n")
        doc_end += 1

    insert_after = doc_end

# Insert Option C block if missing
tag = "SCR->CDM energy-transfer consistency (Option C)"
if any(tag in ln for ln in lines):
    p.write_text("".join(lines))
    print("[OK] Doc block fixed; Option C already present")
    sys.exit(0)

optionc = r'''
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

lines.insert(insert_after, optionc)
p.write_text("".join(lines))
print("[OK] Doc block fixed + Option C inserted")
