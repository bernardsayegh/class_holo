cd /workspaces/class_holo

cat > cpl_class_verify.py << 'PY'
#!/usr/bin/env python3
"""
Run the CPL surrogate parameters through CLASS to explicitly verify
that the background-matched expansion history does NOT reproduce
the holographic S8 suppression.
"""
import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import brentq
import numpy as np

THETA_S_TARGET = 1.040423
c_km_s = 299792.458

BASE_PARAMS = {
    "output": "mPk,tCl",
    "P_k_max_1/Mpc": 10.0,
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

def run_model(label, extra_params, match_theta=True):
    params = dict(BASE_PARAMS)
    params.update(extra_params)
    if match_theta:
        def residual(h):
            p = dict(params)
            p["h"] = float(h)
            c = Class()
            c.set(p)
            c.compute()
            theta = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
            c.struct_cleanup()
            c.empty()
            return theta - THETA_S_TARGET
        h = brentq(residual, 0.55, 0.80, xtol=1e-8, rtol=1e-12)
        params["h"] = float(h)
    c = Class()
    c.set(params)
    c.compute()
    H0 = c.Hubble(0) * c_km_s
    Om = c.Omega_m()
    s8 = c.sigma8()
    S8 = s8 * np.sqrt(Om / 0.3)
    print(f"\n  {label}:")
    print(f"    h       = {params['h']:.5f}")
    print(f"    H0      = {H0:.3f} km/s/Mpc")
    print(f"    Omega_m = {Om:.4f}")
    print(f"    sigma8  = {s8:.4f}")
    print(f"    S8      = {S8:.4f}")
    c.struct_cleanup()
    c.empty()
    return {"H0": H0, "Om": Om, "s8": s8, "S8": S8}

print("=" * 70)
print("CPL SURROGATE: EXPLICIT CLASS VERIFICATION")
print("=" * 70)

lcdm = run_model("LCDM (w=-1)", {})

holo = run_model("Holographic (beta=1/12)", {
    "interaction_beta": 1.0/12.0,
    "f_clust": 0.0,
})

cpl = run_model("CPL surrogate (w0=-1.003, wa=-0.008)", {
    "Omega_fld": 1e-10,
    "use_ppf": "yes",
    "w0_fld": -1.003,
    "wa_fld": -0.008,
})

print("\n" + "=" * 70)
print("COMPARISON")
print("=" * 70)
print(f"  {'Model':<35s} {'sigma8':>8s} {'S8':>8s}")
print(f"  {'-'*35} {'-'*8} {'-'*8}")
print(f"  {'LCDM':<35s} {lcdm['s8']:8.4f} {lcdm['S8']:8.4f}")
print(f"  {'CPL surrogate':<35s} {cpl['s8']:8.4f} {cpl['S8']:8.4f}")
print(f"  {'Holographic (beta=1/12)':<35s} {holo['s8']:8.4f} {holo['S8']:8.4f}")
print(f"\n  CPL vs LCDM:         delta_S8 = {cpl['S8'] - lcdm['S8']:+.4f}")
print(f"  Holographic vs LCDM: delta_S8 = {holo['S8'] - lcdm['S8']:+.4f}")
print("=" * 70)
PY

python3 cpl_class_verify.py
