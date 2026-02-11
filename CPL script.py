cd /workspaces/class_holo

cat > cpl_base.py << 'PY'
#!/usr/bin/env python3
import sys
sys.path.insert(0, "python")
from classy import Class
from scipy.optimize import minimize, brentq
from scipy.integrate import quad
import numpy as np

# -----------------------------
# Targets / controls (consistent with cosmo_stats_thetas_matched.py)
# -----------------------------
THETA_S_TARGET = 1.040423
H_BRACKET = (0.55, 0.80)
BETA_MODEL_A = 1.0/12.0
F_CLUST = 0.0

c_km_s = 299792.458
OMEGA_R = 8.85e-5
Z_FIT = np.linspace(0.0, 2.5, 100)

# Explicit late/early-universe defaults for reproducibility
BASE_PARAMS = {
    "output": "mPk",
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

def get_class_instance(h, beta):
    params = dict(BASE_PARAMS)
    params["h"] = float(h)
    if beta > 0:
        params["interaction_beta"] = float(beta)
        params["f_clust"] = float(F_CLUST)
    c = Class()
    c.set(params)
    c.compute()
    return c

def solve_h_for_theta(beta, bracket=H_BRACKET):
    def f(h):
        c = get_class_instance(h, beta)
        theta = c.get_current_derived_parameters(["100*theta_s"])["100*theta_s"]
        c.struct_cleanup()
        c.empty()
        return theta - THETA_S_TARGET
    return brentq(f, bracket[0], bracket[1], xtol=1e-8, rtol=1e-12)

def E_CPL(z, Om, w0, wa):
    Ode = 1.0 - Om - OMEGA_R
    if Ode <= 0:
        return 1e10
    return np.sqrt(Om*(1+z)**3 + OMEGA_R*(1+z)**4 + Ode*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z)))

def D_M_CPL_arr(z_arr, H0, Om, w0, wa):
    D_M = np.zeros_like(z_arr, dtype=float)
    for i, zi in enumerate(z_arr):
        if zi < 0.001:
            D_M[i] = 0.0
        else:
            result, _ = quad(lambda zp: 1.0/E_CPL(zp, Om, w0, wa), 0, zi, limit=100)
            D_M[i] = c_km_s / H0 * result
    return D_M

def chi2_CPL(params, H_holo, D_M_holo, z_arr):
    H0, Om, w0, wa = params
    if Om <= 0.1 or Om >= 0.9 or H0 <= 50 or H0 >= 90:
        return 1e30
    if 1 - Om - OMEGA_R <= 0:
        return 1e30
    H_cpl = H0 * np.array([E_CPL(z, Om, w0, wa) for z in z_arr])
    D_M_cpl = D_M_CPL_arr(z_arr, H0, Om, w0, wa)
    mask = z_arr > 0.01
    return np.sum(np.log(H_cpl/H_holo)**2) + np.sum(np.log(D_M_cpl[mask]/D_M_holo[mask])**2)

def fit_CPL(H_holo, D_M_holo, z_arr, H0_init):
    best_chi2, best_x = np.inf, None
    for Om0 in [0.26, 0.30, 0.34]:
        for w00 in [-1.05, -1.0, -0.95]:
            for wa0 in [-0.1, 0.0, 0.1]:
                result = minimize(chi2_CPL, [H0_init, Om0, w00, wa0], args=(H_holo, D_M_holo, z_arr),
                                method='Nelder-Mead', options={'maxiter': 3000})
                if result.fun < best_chi2:
                    best_chi2, best_x = result.fun, result.x
    return best_x

def extract_Hz_DM(cosmo, z_arr):
    H = np.array([cosmo.Hubble(z) * c_km_s for z in z_arr])
    D_M = np.array([(1+z) * cosmo.angular_distance(z) for z in z_arr])
    return H, D_M

def compute_residuals(H_holo, D_M_holo, H0, Om, w0, wa, z_arr):
    H_cpl = H0 * np.array([E_CPL(z, Om, w0, wa) for z in z_arr])
    D_M_cpl = D_M_CPL_arr(z_arr, H0, Om, w0, wa)
    mask = z_arr > 0.01
    dH = np.max(np.abs(H_cpl - H_holo) / H_holo) * 100
    dD = np.max(np.abs(D_M_cpl[mask] - D_M_holo[mask]) / D_M_holo[mask]) * 100
    return dH, dD

print("=" * 70)
print("CPL SURROGATE FOR HOLOGRAPHIC MODEL (β=1/12)")
print("=" * 70)

h = solve_h_for_theta(BETA_MODEL_A)
cosmo = get_class_instance(h, BETA_MODEL_A)

H0 = cosmo.Hubble(0) * c_km_s
Om = cosmo.Omega_m()
s8 = cosmo.sigma8()
S8 = s8 * np.sqrt(Om/0.3)

H_holo, D_M_holo = extract_Hz_DM(cosmo, Z_FIT)

x = fit_CPL(H_holo, D_M_holo, Z_FIT, H0)
H0_cpl, Om_cpl, w0, wa = x
dH, dD = compute_residuals(H_holo, D_M_holo, H0_cpl, Om_cpl, w0, wa, Z_FIT)

print(f"\nHolographic Model (β=1/12):")
print(f"  h       = {h:.5f}")
print(f"  H0      = {H0:.3f} km/s/Mpc")
print(f"  Ω_m     = {Om:.3f}")
print(f"  σ8      = {s8:.3f}")
print(f"  S8      = {S8:.3f}")

print(f"\nCPL Surrogate (best-fit):")
print(f"  Ω_m,eff = {Om_cpl:.3f}")
print(f"  w0      = {w0:.3f}")
print(f"  wa      = {wa:.3f}")

print(f"\nFit Quality:")
print(f"  max|δH/H|   = {dH:.4f}%")
print(f"  max|δD_M/D_M| = {dD:.4f}%")

print(f"\nDiscrepancy: Ω_m,true = {Om:.3f}, Ω_m,CPL = {Om_cpl:.3f} ({100*(Om-Om_cpl)/Om:+.1f}%)")
print("=" * 70)

cosmo.struct_cleanup()
cosmo.empty()
PY

python3 cpl_base.py
