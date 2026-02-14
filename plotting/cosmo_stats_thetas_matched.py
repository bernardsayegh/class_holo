import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np

# -----------------------------
# Targets / controls
# -----------------------------
THETA_S_TARGET = 1.040423   # target is 100*theta_s
H_BRACKET = (0.55, 0.80)
BETA_LCDM = 0.0
BETA_MODEL_A = 1.0/12.0
BETA_MODEL_B = 1.0/12.0
F_CLUST = 0.0
PKMAX = 10.0
# -----------------------------

BASE_PARAMS = {
    "output": "mPk",
    "P_k_max_1/Mpc": PKMAX,

    "omega_b": 0.02242,
    "omega_cdm": 0.11933,
    "n_s": 0.9665,
    "A_s": 2.1e-9,
    "tau_reio": 0.054,

    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,

    "Omega_k": 0.0,
}

def compute_stats(h: float, beta: float, super_schw: bool = False, Amap: float = 0.0):
    p = dict(BASE_PARAMS)
    p["h"] = float(h)
    if beta > 0:
        p["interaction_beta"] = float(beta)
        p["f_clust"] = float(F_CLUST)
    if super_schw:
        p["super_schwarzschild_correction"] = "yes"
        p["super_schw_Amap"] = float(Amap)
        p["super_schw_deltaS"] = 0.03
        p["super_schw_gamma"] = 2.0

    c = Class()
    c.set(p)
    c.compute()

    derived = c.get_current_derived_parameters(["100*theta_s"])
    theta = derived["100*theta_s"]
    H0 = 100.0 * c.h()
    Om = float(c.Omega_m())
    sig8 = float(c.sigma8())
    S8 = sig8 * np.sqrt(Om / 0.3)

    # Try to read H0_local (only exists for super-Schwarzschild models)
    H0_local = None
    try:
        ba = c.get_background()
        if "H0_local" in ba:
            H0_local = ba["H0_local"][-1]
        else:
            d2 = c.get_current_derived_parameters(["H0_local"])
            H0_local = d2.get("H0_local", None)
    except:
        pass

    c.struct_cleanup()
    c.empty()
    return float(theta), float(H0), Om, sig8, float(S8), H0_local

def solve_h_for_theta(beta: float, target: float, bracket=(0.55, 0.80),
                      super_schw=False, Amap=0.0):
    a, b = bracket

    def f(h):
        theta, *_ = compute_stats(h, beta, super_schw, Amap)
        return theta - target

    fa = f(a)
    fb = f(b)
    print(f"Bracket check (beta={beta:.6f}): h={a:.4f} -> Δθ={fa:+.6e}, h={b:.4f} -> Δθ={fb:+.6e}")
    if fa * fb > 0:
        raise RuntimeError(
            f"Root not bracketed. Got fa={fa:+.6e}, fb={fb:+.6e}."
        )
    return brentq(f, a, b, xtol=1e-8, rtol=1e-12, maxiter=200)

def run_model(name: str, beta: float, super_schw: bool = False, Amap: float = 0.0):
    h_sol = solve_h_for_theta(beta, THETA_S_TARGET, H_BRACKET, super_schw, Amap)
    theta, H0, Om, sig8, S8, H0_local = compute_stats(h_sol, beta, super_schw, Amap)
    print(f"\n{name}")
    print(f"  beta         = {beta:.12f}")
    print(f"  h (solved)   = {h_sol:.8f}")
    print(f"  100*theta_s  = {theta:.6f}   (target {THETA_S_TARGET:.6f})")
    print(f"  H0           = {H0:.4f} km/s/Mpc")
    if H0_local is not None:
        print(f"  H0_local     = {H0_local:.4f} km/s/Mpc")
    print(f"  Omega_m      = {Om:.8f}")
    print(f"  sigma8       = {sig8:.8f}")
    print(f"  S8           = {S8:.8f}")
    return h_sol, theta, H0, Om, sig8, S8, H0_local

if __name__ == "__main__":
    print("CLASS stats (θs-matched)")
    print("=" * 78)
    print(f"Target: 100*theta_s = {THETA_S_TARGET:.6f}")
    print(f"Bracket: h in [{H_BRACKET[0]}, {H_BRACKET[1]}]")
    print("=" * 78)

    run_model("ΛCDM (β=0)", BETA_LCDM)
    run_model("Model A (β=1/12)", BETA_MODEL_A)
    run_model("Model B (β=1/12, Amap=2)", BETA_MODEL_B, super_schw=True, Amap=2.0)
