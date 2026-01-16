#!/usr/bin/env python3
import re, sys, pathlib
ROOT = pathlib.Path(".").resolve()
bg_h = ROOT/"include/background.h"
bg_c = ROOT/"source/background.c"
pt_c = ROOT/"source/perturbations.c"
for p in (bg_h,bg_c,pt_c):
    if not p.exists(): print("[FATAL] missing",p); sys.exit(1)

def add_struct_indices(t):
    need_rho = "index_bg_rho_scr" not in t
    need_Q   = "index_bg_Q_scr_to_cdm_over_H" not in t
    if not (need_rho or need_Q): return t, False
    m = re.search(r"struct background\s*{", t)
    if not m: raise RuntimeError("no struct background {")
    ins = "\n  /* --- Option C: scr -> cdm transfer (export Q/H to perturbations) --- */\n"
    if need_rho: ins += "  int index_bg_rho_scr;\n"
    if need_Q:   ins += "  int index_bg_Q_scr_to_cdm_over_H;\n"
    pos = m.end()
    return t[:pos]+ins+t[pos:], True

def alloc_indices(t):
    # anchor on rho_cdm allocation if present, else rho_tot
    m = re.search(r"(pba->index_bg_rho_cdm\s*=\s*index_bg;\s*index_bg\+\+;)", t) \
        or re.search(r"(pba->index_bg_rho_tot\s*=\s*index_bg;\s*index_bg\+\+;)", t)
    if not m: raise RuntimeError("cannot find index allocation anchor (rho_cdm or rho_tot)")
    ins = ""
    if "index_bg_rho_scr" not in t: ins += "\n  pba->index_bg_rho_scr = index_bg; index_bg++;"
    if "index_bg_Q_scr_to_cdm_over_H" not in t: ins += "\n  pba->index_bg_Q_scr_to_cdm_over_H = index_bg; index_bg++;"
    if not ins: return t, False
    return t.replace(m.group(1), m.group(1) + "\n  /* --- Option C: scr -> cdm transfer --- */" + ins + "\n", 1), True

def export_pvecback(t):
    if "pvecback[pba->index_bg_Q_scr_to_cdm_over_H]" in t: return t, False
    am = re.search(r"(pvecback\[pba->index_bg_rho_cdm\]\s*=\s*[^;]+;)", t) \
        or re.search(r"(pvecback\[pba->index_bg_rho_tot\]\s*=\s*[^;]+;)", t)
    if not am: raise RuntimeError("cannot find pvecback anchor (rho_cdm or rho_tot)")
    if "rho_scr" not in t or "Q_scr_to_cdm_over_H" not in t:
        raise RuntimeError("background.c missing symbols rho_scr and/or Q_scr_to_cdm_over_H (rename patcher if yours differ)")
    ins = ("\n  /* --- Option C: export scr reservoir + Q/H to perturbations --- */\n"
           "  pvecback[pba->index_bg_rho_scr] = rho_scr;\n"
           "  pvecback[pba->index_bg_Q_scr_to_cdm_over_H] = Q_scr_to_cdm_over_H;\n")
    return t.replace(am.group(1), am.group(1)+ins, 1), True

def patch_pert(t):
    if "SCR->CDM energy-transfer consistency (Option C)" in t: return t, False
    if "BEGIN HOLOGRAPHIC INTERACTION" not in t:
        raise RuntimeError("no BEGIN HOLOGRAPHIC INTERACTION marker in perturbations.c")
    block = r"""
/* --- SCR->CDM energy-transfer consistency (Option C) -------------------------
   Use the SAME Q_scr_to_cdm_over_H(a) as background.c.
   Frame: Q^μ || u_cdm (no momentum transfer in CDM frame).
   Default δQ choice: δQ/Q = δ_cdm.
------------------------------------------------------------------------------- */
if (pba->super_schw_gamma != 0.) {
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double rho_cdm = pvecback[pba->index_bg_rho_cdm];
  double Qscr_over_H = pvecback[pba->index_bg_Q_scr_to_cdm_over_H];

  if (rho_cdm > 0. && Qscr_over_H != 0.) {
    double coupling = a * H * Qscr_over_H / rho_cdm;
    double delta_cdm = y[pv->index_pt_delta_cdm];
    double deltaQ_over_Q = delta_cdm; /* δQ/Q = δ_cdm */
    dy[pv->index_pt_delta_cdm] += coupling * (delta_cdm - deltaQ_over_Q);
  }
}
"""
    return re.sub(r"(\/\*\s*BEGIN HOLOGRAPHIC INTERACTION\s*\*\/)", r"\1"+block, t, count=1), True

h = bg_h.read_text(); h2, ch = add_struct_indices(h)
if ch: bg_h.write_text(h2); print("[OK] background.h indices ensured")
else: print("[SKIP] background.h ok")

c = bg_c.read_text(); c2, ch = alloc_indices(c)
if ch: bg_c.write_text(c2); print("[OK] background.c alloc ensured")
else: print("[SKIP] background.c alloc ok")

c = bg_c.read_text(); c3, ch = export_pvecback(c)
if ch: bg_c.write_text(c3); print("[OK] background.c pvecback export ensured")
else: print("[SKIP] background.c pvecback export ok")

p = pt_c.read_text(); p2, ch = patch_pert(p)
if ch: pt_c.write_text(p2); print("[OK] perturbations.c patched")
else: print("[SKIP] perturbations.c already patched")
