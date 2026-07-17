#!/usr/bin/env python3
import os, glob, re, math
import numpy as np
BURN_FRAC = 0.30
ROOT = os.path.expanduser("~/class_holo_test")
os.chdir(ROOT)
# ---------- v4 TARGETS ----------
H0_TARGET, H0_TARGET_ERR = 73.49, 0.93        # Riess et al. 2025 (ApJL 992, L34), HST+JWST Cepheids
S8_TARGET, S8_TARGET_ERR = 0.776, 0.017
# ---------- RUNS ----------
# (path, label, Δk vs ΛCDM, per-run H0 target or None -> default R25)
runs = [
    ("chains_test/lcdm_v5_acc",              "ΛCDM",                0, None),
    ("chains_test/lcdm_v5_acc_noDES",        "ΛCDM (noDES)",        0, None),
    ("chains_test/lcdm_noprior_v5_acc",      "ΛCDM (noH0)",         0, None),
    ("chains_test/lcdm_A2_nd_v5_acc",        "ΛCDM, A2, ND",        0, None),
    ("chains_test/lcdm_Afree_nd_v5_acc",     "ΛCDM, Af, ND",        1, None),
    ("chains_test/modelA_v5_acc",            "Model A (β=1/12)",    0, None),
    ("chains_test/modelA_bfree_v5_acc",      "Model A (β free)",    1, None),
    ("chains_test/modelA_fc1_v5_acc",        "Model A (fc=1)",      0, None),
    ("chains_test/modelA_fcfree_v5_acc",     "Model A (fc free)",   1, None),
    ("chains_test/modelB_A1_v5_acc",         "Model B (A=1)",       0, None),
    ("chains_test/modelB_A2_v5_acc",         "Model B (A=2)",       0, None),
    ("chains_test/modelB_A2_v5_acc_noDES",   "Model B (A2,noDES)",  0, None),
    ("chains_test/modelB_Afree_v5_acc",      "Model B (Af)",        1, None),
    ("chains_test/modelB_A2_noprior_v5_acc", "Model B (A2,noH0)",   0, None),
    ("chains_test/modelB_A2_trgb_v5_acc",    "Model B (A2,TRGB)",   0, (69.96, 1.53)),
    ("chains_test/modelAB_bfree_Afree_v5_acc", "Model AB (β,A free)", 2, None),
    ("chains_test/modelD_ampfree_v5_acc",     "Model D (amp free)",  1, None),
    ("chains_test/modelC_ampfree_Afree_v5_acc","Model C (amp,Af)",    2, None),
    ("chains_test/modelC_v5_acc",            "Model C (resv+map)",  0, None),
    # v3 taxonomy row: fit under SH0ES 73.04±1.04 -> Δχ² NOT comparable to v4 ΛCDM.
    # Re-run under R25 before adding here:
    # ("chains_test/modelB_A2_ode1_v3_full", "Model B (A2,ODE1)[v3]", 0, None),
]
LCDM_NAME = "ΛCDM"
# ---------- UTILITIES ----------
def has_chain(run):
    return bool(glob.glob(f"{run}.*.txt"))
def load_raw_chain(run, burn_frac=BURN_FRAC):
    fns = sorted(glob.glob(f"{run}.*.txt"))
    fns = [f for f in fns if f.split(".")[-2].isdigit()]
    if not fns:
        return None
    with open(fns[0]) as f:
        hdr = f.readline().lstrip("#").split()
    names = hdr[2:]
    chunks = []
    for fn in fns:
        try:
            chunks.append(np.loadtxt(fn))
        except Exception:
            pass
    if not chunks:
        return None
    chunks = [c[int(burn_frac*len(c)):] for c in chunks if c.ndim == 2 and len(c) >= 10]
    if not chunks:
        return None
    data = np.vstack(chunks)
    if len(data) < 10:
        return None
    ncols = data.shape[1]
    offset = ncols - len(names)
    if offset < 2:
        raise RuntimeError(f"Unexpected column layout: {ncols} cols, {len(names)} params")
    colmap = {n: i + offset for i, n in enumerate(names)}
    raw = {"run": run, "fn": fns[0], "names": names, "data": data, "col": colmap}
    for k in ("H0_local", "S8", "S8_loc"):
        if k not in raw["col"]:
            return raw
    H = raw["data"][:, raw["col"]["H0_local"]]
    S = raw["data"][:, raw["col"]["S8"]]
    Sl = raw["data"][:, raw["col"]["S8_loc"]]
    mH, mS, mSl = float(np.mean(H)), float(np.mean(S)), float(np.mean(Sl))
    if (mH < 5.0) and (mS > 5.0) and (mSl < 5.0):
        raw["col"]["H0_local"], raw["col"]["S8"], raw["col"]["S8_loc"] = (
            raw["col"]["S8"], raw["col"]["S8_loc"], raw["col"]["H0_local"]
        )
    return raw
def col(raw, name):
    j = raw["col"].get(name, None)
    return raw["data"][:, j] if j is not None else None
def wmean_wstd(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    s = np.sum(w)
    if s <= 0:
        return (None, None)
    mu = float(np.sum(w * x) / s)
    var = float(np.sum(w * (x - mu)**2) / s)
    return (mu, float(np.sqrt(max(var, 0.0))))
def bestfit_idx(raw):
    return int(np.argmin(raw["data"][:, 1]))
def global_S8(raw):
    w = raw["data"][:, 0]
    S8 = col(raw, "S8")
    if S8 is not None:
        return wmean_wstd(S8, w)
    sig8 = col(raw, "sigma8"); Om = col(raw, "Omega_m")
    if sig8 is None or Om is None:
        return (None, None)
    return wmean_wstd(sig8 * np.sqrt(Om / 0.3), w)
def pick_h0loc_s8loc(raw):
    w = raw["data"][:, 0]
    H0_local = col(raw, "H0_local"); S8_loc = col(raw, "S8_loc")
    h0m = h0s = s8m = s8s = None
    if H0_local is not None:
        h0m, h0s = wmean_wstd(H0_local, w)
    if S8_loc is not None:
        s8m, s8s = wmean_wstd(S8_loc, w)
    return (h0m, h0s, s8m, s8s)
# ---------- CHI2 ----------
def chi2_total_bestfit(raw):
    i = bestfit_idx(raw)
    for k in ["chi2", "chi2__total", "chi2__chi2"]:
        a = col(raw, k)
        if a is not None:
            return (float(a[i]), True)
    tot = 0.0; ok = False
    seen_highl = seen_lowee = seen_lens = False
    for k in [
        "chi2__planck_NPIPE_highl_CamSpec.TTTEEE",
        "chi2__planck_2018_highl_plik.TTTEEE",
        "chi2__planck_2018_lowl.TT",
        "chi2__planck_2018_lowl.EE_sroll2",
        "chi2__planck_2018_lowl.EE",
        "chi2__planck_2018_lensing.CMBMarged",
        "chi2__planck_2018_lensing.clik",
    ]:
        a = col(raw, k)
        if a is not None:
            if "highl" in k and seen_highl: continue
            if "lowl.EE" in k and seen_lowee: continue
            if "lensing" in k and seen_lens: continue
            tot += float(a[i]); ok = True
            if "highl" in k: seen_highl = True
            if "lowl.EE" in k: seen_lowee = True
            if "lensing" in k: seen_lens = True
    for k in ["chi2__bao.desi_dr2",
              "chi2__bao.sixdf_2011_bao", "chi2__bao.sdss_dr7_mgs",
              "chi2__bao.sdss_dr12_consensus_full_shape",
              "chi2__bao.sdss_dr16_baoplus_elg", "chi2__bao.sdss_dr16_baoplus_qso",
              "chi2__bao.sdss_dr16_baoplus_lyauto", "chi2__bao.sdss_dr16_baoplus_lyxqso"]:
        a = col(raw, k)
        if a is not None:
            tot += float(a[i]); ok = True
    a = col(raw, "chi2__sn.pantheonplus")
    if a is not None:
        tot += float(a[i]); ok = True
    for k in ["chi2__desy3_likelihood.DESY3",
              "chi2__desy3_likelihood_fast.DESY3",
              "chi2__desy3_likelihood_h0local.DESY3",
              "chi2__desy3_likelihood_scrloc.DESY3"]:
        a = col(raw, k)
        if a is not None:
            tot += float(a[i]); ok = True; break
    for k in ["chi2__scrboost.SCRBoost_SH0ES",
              "chi2__shoes_h0local.SH0ES_H0local", "chi2__shoes_h0local.SH0ES_H0",
              "chi2__shoes_h0.SH0ES_H0", "chi2__trgb_h0local.TRGB_H0local"]:
        a = col(raw, k)
        if a is not None:
            tot += float(a[i]); ok = True; break
    return (tot if ok else float("nan"), ok)
def chi2_shoes_bestfit(raw):
    i = bestfit_idx(raw)
    for k in ["chi2__scrboost.SCRBoost_SH0ES",
              "chi2__shoes_h0local.SH0ES_H0local", "chi2__shoes_h0local.SH0ES_H0",
              "chi2__shoes_h0.SH0ES_H0", "chi2__trgb_h0local.TRGB_H0local"]:
        a = col(raw, k)
        if a is not None:
            return float(a[i])
    return 0.0
# ---------- R-1 ----------
_R_EXPLICIT = [
    re.compile(r"\bR\s*-\s*1\b\s*[:=]\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)"),
    re.compile(r"\bRminus1\b\s*[:=]\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", re.IGNORECASE),
]
def read_rminus1_from_progress(run):
    prog = f"{run}.progress"
    if not os.path.exists(prog):
        return None
    try:
        lines = [ln.rstrip("\n") for ln in open(prog, "r", errors="ignore") if ln.strip()]
        for ln in reversed(lines[-2000:]):
            for pat in _R_EXPLICIT:
                m = pat.search(ln)
                if m:
                    return m.group(1)
        for ln in reversed(lines[-2000:]):
            parts = ln.split()
            if len(parts) >= 4:
                try:
                    r1 = float(parts[3])
                    if np.isfinite(r1):
                        return str(r1)
                except Exception:
                    pass
    except Exception:
        pass
    return None
def split_rminus1_estimate(raw):
    keys = ["H0", "H0_local", "Omega_m", "sigma8", "S8", "n_s", "tau_reio"]
    present = [k for k in keys if k in raw["col"]]
    if not present:
        return None
    N = raw["data"].shape[0]
    if N < 80:
        return None
    n = N // 2
    rminus = []
    for k in present:
        a = col(raw, k)
        if a is None:
            continue
        x1 = a[:n]; x2 = a[n:2*n]
        mean1, mean2 = float(np.mean(x1)), float(np.mean(x2))
        mean = 0.5 * (mean1 + mean2)
        var1, var2 = float(np.var(x1, ddof=1)), float(np.var(x2, ddof=1))
        W = 0.5 * (var1 + var2)
        if W <= 0:
            continue
        B = n * ((mean1 - mean)**2 + (mean2 - mean)**2)
        Var_hat = (n - 1.0) / n * W + (1.0 / n) * B
        if Var_hat <= 0:
            continue
        rm = math.sqrt(Var_hat / W) - 1.0
        if np.isfinite(rm):
            rminus.append(rm)
    return float(max(0.0, np.max(rminus))) if rminus else None
# ---------- PER-MODEL BLOCK ----------
def print_model_block(raw, name):
    w = raw["data"][:, 0]
    i = bestfit_idx(raw)
    def show(label, key, nd=6):
        a = col(raw, key)
        if a is None:
            return
        mu, sd = wmean_wstd(a, w)
        print(f"  {label:<18}: {mu:.{nd}f} ± {sd:.{nd}f}")
    print("\n" + "="*90)
    print(f"{name} (N={raw['data'].shape[0]})")
    print("="*90)
    print("\nParameters (mean ± std):")
    for label, key in [
        ("H0", "H0"), ("H0_local", "H0_local"),
        ("omega_b", "omega_b"), ("omega_cdm", "omega_cdm"),
        ("sigma8", "sigma8"), ("Omega_m", "Omega_m"),
        ("n_s", "n_s"), ("tau_reio", "tau_reio"),
        ("interaction_beta", "interaction_beta"),
        ("super_schw_Amap", "super_schw_Amap"),
        ("super_schw_amp", "super_schw_amp"),
        ("f_clust", "f_clust"),
        ("A_IA", "A_IA"), ("alpha_IA", "alpha_IA"),
    ]:
        show(label, key)
    s8m, s8s = global_S8(raw)
    if s8m is not None:
        print(f"  {'S8 (derived)':<18}: {s8m:.6f} ± {s8s:.6f}")
    h0m, h0s, s8lm, s8ls = pick_h0loc_s8loc(raw)
    if h0m is not None:
        print(f"  {'H0_loc':<18}: {h0m:.6f} ± {h0s:.6f}")
    if s8lm is not None:
        print(f"  {'S8_loc':<18}: {s8lm:.6f} ± {s8ls:.6f}")
    print("\nBest-fit χ² breakdown:")
    def bf(key):
        a = col(raw, key)
        return float(a[i]) if a is not None else None
    chi2_cmb = 0.0
    for k, lab in [
        ("chi2__planck_NPIPE_highl_CamSpec.TTTEEE", "NPIPE CamSpec TTTEEE"),
        ("chi2__planck_2018_highl_plik.TTTEEE", "highl_plik.TTTEEE"),
        ("chi2__planck_2018_lowl.TT", "lowl.TT"),
        ("chi2__planck_2018_lowl.EE_sroll2", "lowl.EE_sroll2"),
        ("chi2__planck_2018_lowl.EE", "lowl.EE"),
        ("chi2__planck_2018_lensing.CMBMarged", "lensing.CMBMarged"),
        ("chi2__planck_2018_lensing.clik", "lensing.clik"),
    ]:
        v = bf(k)
        if v is not None:
            chi2_cmb += v
            print(f"  {lab:<35}: {v:.1f}")
    if chi2_cmb > 0:
        print(f"  {'CMB TOTAL':<35}: {chi2_cmb:.1f}")
    bao = 0.0; any_bao = False
    v = bf("chi2__bao.desi_dr2")
    if v is not None:
        bao += v; any_bao = True
        print(f"  {'DESI DR2':<35}: {v:.2f}")
    for k, lab in [
        ("chi2__bao.sixdf_2011_bao", "  6dFGS"),
        ("chi2__bao.sdss_dr7_mgs", "  DR7 MGS"),
        ("chi2__bao.sdss_dr12_consensus_full_shape", "  DR12 consensus"),
        ("chi2__bao.sdss_dr16_baoplus_elg", "  eBOSS ELG"),
        ("chi2__bao.sdss_dr16_baoplus_qso", "  eBOSS QSO"),
        ("chi2__bao.sdss_dr16_baoplus_lyauto", "  Ly-a auto"),
        ("chi2__bao.sdss_dr16_baoplus_lyxqso", "  Ly-a x QSO"),
    ]:
        v = bf(k)
        if v is not None:
            bao += v; any_bao = True
            print(f"  {lab:<35}: {v:.2f}")
    if any_bao:
        print(f"  {'BAO TOTAL':<35}: {bao:.1f}")
    v = bf("chi2__sn.pantheonplus")
    if v is not None:
        print(f"  {'Pantheon+':<35}: {v:.1f}")
    for k in ["chi2__desy3_likelihood.DESY3",
              "chi2__desy3_likelihood_fast.DESY3",
              "chi2__desy3_likelihood_h0local.DESY3",
              "chi2__desy3_likelihood_scrloc.DESY3"]:
        v = bf(k)
        if v is not None:
            print(f"  {'DES Y3':<35}: {v:.1f}"); break
    for k, lab in [
        ("chi2__scrboost.SCRBoost_SH0ES", "SH0ES (SCRBoost)"),
        ("chi2__shoes_h0local.SH0ES_H0local", "SH0ES (H0_local)"),
        ("chi2__shoes_h0local.SH0ES_H0", "SH0ES (H0_local)"),
        ("chi2__shoes_h0.SH0ES_H0", "SH0ES (physical H0)"),
        ("chi2__trgb_h0local.TRGB_H0local", "TRGB (H0_local)"),
    ]:
        v = bf(k)
        if v is not None:
            print(f"  {lab:<35}: {v:.1f}"); break
    tot, ok = chi2_total_bestfit(raw)
    if ok:
        print(f"  {'TOTAL':<35}: {tot:.1f}")
    print("\nH0 summary:")
    H0 = col(raw, "H0")
    if H0 is not None:
        H0m, H0s = wmean_wstd(H0, w)
        print(f"  H0 (physical): {H0m:6.2f} ± {H0s:.2f}")
    h0m, h0s, _, _ = pick_h0loc_s8loc(raw)
    if h0m is not None:
        print(f"  H0_loc:        {h0m:6.2f} ± {h0s:.2f}")
    r1 = read_rminus1_from_progress(raw["run"])
    if r1 is None:
        est = split_rminus1_estimate(raw)
        r1 = f"{est:.6f}*" if est is not None else "—"
    print(f"  R-1:           {r1}")
# ========== MAIN ==========
results = []
lcdm_chi2 = None
for run, name, dk, tgt in runs:
    if not has_chain(run):
        continue
    raw = load_raw_chain(run, BURN_FRAC)
    if raw is None:
        continue
    print_model_block(raw, name)
    w = raw["data"][:, 0]
    n = raw["data"].shape[0]
    H0 = col(raw, "H0")
    H0m, H0s = (wmean_wstd(H0, w) if H0 is not None else (None, None))
    h0loc_m, h0loc_s, s8loc_m, s8loc_s = pick_h0loc_s8loc(raw)
    s8m, s8s = global_S8(raw)
    chi2, ok = chi2_total_bestfit(raw)
    shoes = chi2_shoes_bestfit(raw)
    r1 = read_rminus1_from_progress(run)
    if r1 is None:
        est = split_rminus1_estimate(raw)
        r1 = f"{est:.6f}*" if est is not None else "—"
    results.append({
        "name": name, "n": n, "chi2": chi2, "shoes": shoes, "dk": dk,
        "H0": H0m, "H0s": H0s,
        "H0loc": h0loc_m, "H0locs": h0loc_s,
        "S8": s8m, "S8s": s8s,
        "S8loc": s8loc_m, "S8locs": s8loc_s,
        "r1": r1, "tgt": tgt
    })
    if name == LCDM_NAME:
        lcdm_chi2 = chi2
# ========== SUMMARY TABLE ==========
if results:
    W = 160
    print("\n" + "=" * W)
    print("SUMMARY TABLE v5 (accumulator, Planck nu)  [Planck + DESI DR2 + Pantheon+ + DES Y3 + SH0ES R25]")
    print("=" * W)
    print(f"{'Model':<22} {'N':>6} {'Δk':>3} {'χ²_bf':>8} {'Δχ²':>7} {'ΔAIC':>7} {'SH0ES':>6}"
          f" {'H0_phys':>12} {'H0_loc':>12} {'S8':>14} {'S8_loc':>14}"
          f" {'H0 σ':>6} {'S8 σ':>6} {'R-1':>12}")
    print("-" * W)
    for r in results:
        def fmt(mu, sd, nd):
            return f"{mu:.{nd}f}±{sd:.{nd}f}" if (mu is not None and sd is not None) else "—"
        H0p = fmt(r["H0"], r["H0s"], 2)
        H0l = fmt(r["H0loc"], r["H0locs"], 2)
        S8s = fmt(r["S8"], r["S8s"], 4)
        S8l = fmt(r["S8loc"], r["S8locs"], 4)
        t_mu, t_err = r["tgt"] if r["tgt"] is not None else (H0_TARGET, H0_TARGET_ERR)
        H0sig = "—"
        if r["H0loc"] is not None:
            H0sig = f"{abs(r['H0loc'] - t_mu) / t_err:.2f}"
        elif r["H0"] is not None:
            H0sig = f"{abs(r['H0'] - t_mu) / t_err:.2f}"
        S8sig = "—"
        if r["S8"] is not None:
            S8sig = f"{abs(r['S8'] - S8_TARGET) / S8_TARGET_ERR:.2f}"
        dchi2_str = daic_str = "—"
        if lcdm_chi2 is not None:
            dchi2 = r["chi2"] - lcdm_chi2
            daic = dchi2 + 2 * r["dk"]
            if r["name"] == LCDM_NAME:
                dchi2_str = "0.0"
                daic_str = "0.0"
            else:
                dchi2_str = f"{dchi2:+.1f}"
                daic_str = f"{daic:+.1f}"
        flag = "*" if r["tgt"] is not None else " "
        print(f"{r['name']:<22} {r['n']:>6} {r['dk']:>3} {r['chi2']:>8.1f} {dchi2_str:>7} {daic_str:>7} {r['shoes']:>6.1f}"
              f" {H0p:>12} {H0l:>12} {S8s:>14} {S8l:>14}"
              f" {H0sig:>5}{flag} {S8sig:>6} {r['r1']:>12}")
    print("-" * W)
    print(f"{'TARGETS (R25)':<22} {'':>6} {'':>3} {'':>8} {'':>7} {'':>7} {'0':>6}"
          f" {'':>12} {'73.49±0.93':>12} {'0.776±0.017':>14} {'0.776±0.017':>14}"
          f" {'0':>6} {'0':>6}")
    print("  * = H0 σ computed against this run's own anchor (TRGB row: CCHP 69.96±1.53)")
else:
    print("No v5 chains found yet — run again once the first checkpoints land.")
