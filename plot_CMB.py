import numpy as np
import matplotlib.pyplot as plt
import glob

# Find the files
lcdm_files = glob.glob('output/lcdm_*_cl_lensed.dat')
holo_files = glob.glob('output/holo_*_cl_lensed.dat')

print(f"LCDM file: {lcdm_files}")
print(f"Holo file: {holo_files}")

lcdm_cl = np.loadtxt(lcdm_files[0])
holo_cl = np.loadtxt(holo_files[0])

# Create figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

ell_lcdm = lcdm_cl[:, 0]
ell_holo = holo_cl[:, 0]
TT_lcdm = lcdm_cl[:, 1]
TT_holo = holo_cl[:, 1]

ax1.semilogy(ell_lcdm, TT_lcdm, 'k-', label=r'$\Lambda$CDM', linewidth=1.5)
ax1.semilogy(ell_holo, TT_holo, 'b--', label=r'Holographic ($\beta=0.5$)', linewidth=1.5)
ax1.set_ylabel(r'$D_\ell^{TT}$ [$\mu$K$^2$]', fontsize=12)
ax1.legend(fontsize=11)
ax1.set_xlim(2, 2500)
ax1.grid(True, alpha=0.3)
ax1.set_title('CMB Temperature Power Spectrum', fontsize=14)

# Ratio plot
min_len = min(len(ell_lcdm), len(ell_holo))
ratio = TT_holo[:min_len] / TT_lcdm[:min_len]

ax2.plot(ell_lcdm[:min_len], ratio, 'b-', linewidth=1.5)
ax2.axhline(1, color='gray', linestyle=':', linewidth=1)
ax2.axhspan(0.99, 1.01, alpha=0.2, color='green', label=r'$\pm 1\%$ band')
ax2.set_xlabel(r'Multipole $\ell$', fontsize=12)
ax2.set_ylabel(r'Holo / $\Lambda$CDM', fontsize=12)
ax2.set_ylim(0.98, 1.02)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

mean_dev = np.mean(np.abs(ratio - 1)) * 100
max_dev = np.max(np.abs(ratio - 1)) * 100
ax2.text(0.02, 0.98, f'Mean deviation: {mean_dev:.3f}%\nMax deviation: {max_dev:.2f}%', 
         transform=ax2.transAxes, fontsize=10, va='top', 
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('Figure_CMB_holographic.png', dpi=150, bbox_inches='tight')
plt.savefig('Figure_CMB_holographic.pdf', bbox_inches='tight')
print(f"CMB comparison: mean dev = {mean_dev:.3f}%, max dev = {max_dev:.2f}%")
print("Saved CMB figures!")
