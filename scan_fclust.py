import sys
sys.path.insert(0, 'python')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

print("Scanning f_clust for two Ω_m scenarios (perturbation only)...")

cases = [
    {"name": "Planck $\\Omega_m$ = 0.31", "omega_cdm": 0.11933, "Omega_m": 0.31, "color": "blue"},
    {"name": "Late Universe $\\Omega_m$ = 1/3", "omega_cdm": 0.1302, "Omega_m": 1/3, "color": "red"},
]

f_clust_values = np.linspace(0.0, 1.0, 21)
results = {}

for case in cases:
    print(f"\n{case['name']}:")
    S8_list = []
    
    for f_clust in f_clust_values:
        cosmo = Class()
        cosmo.set({
            'output': 'mPk',
            'P_k_max_1/Mpc': 10.0,
            'h': 0.6766,
            'omega_b': 0.02242,
            'omega_cdm': case['omega_cdm'],
            'n_s': 0.9665,
            'A_s': 2.1e-9,
            'interaction_beta': 0.5,
            'f_clust': f_clust,
        })
        try:
            cosmo.compute()
            sig8 = cosmo.sigma8()
            S8 = sig8 * np.sqrt(case['Omega_m'] / 0.3)
            S8_list.append(S8)
            print(f"  f_clust={f_clust:.2f}: S8={S8:.4f}")
        except:
            S8_list.append(np.nan)
        cosmo.struct_cleanup()
        cosmo.empty()
    
    results[case['name']] = {'f_clust': f_clust_values, 'S8': np.array(S8_list), 'color': case['color']}

# Find intersections
def find_intersection(f_arr, S8_arr, target):
    for i in range(len(S8_arr)-1):
        if S8_arr[i] >= target >= S8_arr[i+1]:
            f1, f2 = f_arr[i], f_arr[i+1]
            s1, s2 = S8_arr[i], S8_arr[i+1]
            return f1 + (target - s1) * (f2 - f1) / (s2 - s1)
    return None

# Plot
fig, ax = plt.subplots(figsize=(10, 7))

for name, data in results.items():
    ax.plot(data['f_clust'], data['S8'], 'o-', color=data['color'], 
            label=name, linewidth=2, markersize=6)

ax.axhline(y=0.77, color='green', linestyle='--', linewidth=2, label='Weak lensing $S_8$ = 0.77')
ax.axhline(y=0.83, color='orange', linestyle='--', linewidth=2, label='Planck ΛCDM $S_8$ = 0.83')

# Find trajectory points
f_start = find_intersection(results["Planck $\\Omega_m$ = 0.31"]['f_clust'], 
                            results["Planck $\\Omega_m$ = 0.31"]['S8'], 0.83)
f_end = find_intersection(results["Late Universe $\\Omega_m$ = 1/3"]['f_clust'],
                          results["Late Universe $\\Omega_m$ = 1/3"]['S8'], 0.77)

if f_start and f_end:
    ax.annotate('', xy=(f_end, 0.77), xytext=(f_start, 0.83),
                arrowprops=dict(arrowstyle='->', color='black', lw=3))
    ax.plot([f_start], [0.83], 'ko', markersize=10)
    ax.plot([f_end], [0.77], 'ko', markersize=10)
    print(f"\nTrajectory: f_clust {f_start:.3f} -> {f_end:.3f}")

ax.set_xlabel('$f_{\\rm clust}$ (clustering fraction)', fontsize=14)
ax.set_ylabel('$S_8 = \\sigma_8 (\\Omega_m/0.3)^{0.5}$', fontsize=14)
ax.set_title('Holographic Dark Energy: $S_8$ vs Clustering Fraction ($\\beta$ = 0.5)', fontsize=14)
ax.legend(loc='upper right', fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.55, 0.95)

plt.tight_layout()
plt.savefig('S8_fclust.png', dpi=150, bbox_inches='tight')
plt.savefig('S8_fclust.pdf', bbox_inches='tight')
print("\nSaved: S8_fclust.png and S8_fclust.pdf")
