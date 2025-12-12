import numpy as np
import matplotlib.pyplot as plt
import subprocess
import glob
import os

betas = [0.0, 0.25, 0.5, 0.75, 1.0]
results = {}
sigma8_values = {}

# Clean old outputs
for f in glob.glob('output/holo_*'):
    os.remove(f)

for beta in betas:
    root = f'output/holo_b{int(beta*100):03d}_'
    ini = f"""h = 0.6766
omega_b = 0.02242
omega_cdm = 0.11933
n_s = 0.9665
A_s = 2.105e-9
tau_reio = 0.0561
interaction_beta = {beta}
output = mPk
P_k_max_1/Mpc = 10.0
root = {root}
"""
    with open('run_beta.ini', 'w') as f:
        f.write(ini)
    
    # Run CLASS and capture sigma8
    result = subprocess.run(['./class', 'run_beta.ini'], capture_output=True, text=True)
    for line in result.stdout.split('\n'):
        if 'sigma8' in line:
            try:
                sigma8_values[beta] = float(line.split('=')[1].split()[0])
            except:
                pass
    
    # Find the pk file
    pk_files = glob.glob(f'{root}*pk.dat')
    if pk_files:
        data = np.loadtxt(pk_files[0])
        results[beta] = data
        s8 = sigma8_values.get(beta, 0)
        print(f"beta={beta}: sigma8={s8:.4f}, {len(data)} k-points")

# Calculate S8 values
Omega_m = 0.3111
S8_values = {b: s * np.sqrt(Omega_m/0.3) for b, s in sigma8_values.items()}

print("\n=== Summary ===")
for beta in betas:
    if beta in sigma8_values:
        print(f"beta={beta}: sigma8={sigma8_values[beta]:.4f}, S8={S8_values[beta]:.4f}")

# Plot P(k)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

colors = ['black', 'blue', 'red', 'green', 'purple']
labels = [r'$\Lambda$CDM ($\beta=0$)', r'$\beta=0.25$', r'$\beta=0.5$', r'$\beta=0.75$', r'$\beta=1.0$']
styles = ['--', '-', '-', '-', '-']

for i, beta in enumerate(betas):
    if beta in results:
        k = results[beta][:, 0]
        Pk = results[beta][:, 1]
        lbl = labels[i]
        if beta in sigma8_values:
            lbl += f' ($\\sigma_8$={sigma8_values[beta]:.3f})'
        ax1.loglog(k, Pk, styles[i], color=colors[i], label=lbl, linewidth=2)
        
        if beta > 0 and 0.0 in results:
            k_ref = results[0.0][:, 0]
            Pk_ref = results[0.0][:, 1]
            Pk_interp = np.interp(k_ref, k, Pk)
            ratio = Pk_interp / Pk_ref
            ax2.semilogx(k_ref, ratio, styles[i], color=colors[i], linewidth=2)

ax1.set_ylabel(r'$P(k)$ [Mpc$^3$/h$^3$]', fontsize=12)
ax1.legend(fontsize=10, loc='lower left')
ax1.set_xlim(1e-3, 10)
ax1.set_ylim(1, 1e5)
ax1.grid(True, alpha=0.3)
ax1.set_title('Holographic Dark Drag: Matter Power Spectrum', fontsize=14)

ax2.set_xlabel(r'$k$ [h/Mpc]', fontsize=12)
ax2.set_ylabel(r'$P(k)/P_{\Lambda\mathrm{CDM}}$', fontsize=12)
ax2.axhline(1, color='gray', linestyle=':', linewidth=1)
ax2.axhspan(0.88, 0.95, alpha=0.2, color='orange', label='Weak lensing preferred')
ax2.set_ylim(0.70, 1.05)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10, loc='lower left')

plt.tight_layout()
plt.savefig('Figure_Pk_holographic.png', dpi=150, bbox_inches='tight')
plt.savefig('Figure_Pk_holographic.pdf', bbox_inches='tight')
print("\nSaved P(k) figures!")
