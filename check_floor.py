import numpy as np
import subprocess
import glob

# Run CLASS with background output
ini = """h = 0.6766
omega_b = 0.02242
omega_cdm = 0.11933
n_s = 0.9665
A_s = 2.105e-9
tau_reio = 0.0561
interaction_beta = 0.5
output = mPk, background
write background = yes
root = output/floor_test_
"""

with open('floor_test.ini', 'w') as f:
    f.write(ini)

subprocess.run(['./class', 'floor_test.ini'], capture_output=True)

# Find and load background file
bg_files = glob.glob('output/floor_test_*background.dat')
if bg_files:
    print(f"Loading: {bg_files[0]}")
    # Read header to find column names
    with open(bg_files[0], 'r') as f:
        for line in f:
            if line.startswith('#'):
                if 'columns' in line.lower() or ':' in line:
                    print(line.strip())
            else:
                break
    
    # Load data
    data = np.loadtxt(bg_files[0])
    print(f"\nShape: {data.shape}")
    print(f"Columns available: {data.shape[1]}")
    
    # Usually: z, conformal time, H, Omega_m, Omega_r, Omega_Lambda, ...
    # Let's print first few rows
    print("\nFirst 5 rows (first 8 columns):")
    print(data[:5, :8])
    
    print("\nLast 5 rows (z -> 0):")
    print(data[-5:, :8])
else:
    print("No background file found")
    print(glob.glob('output/floor_test_*'))
