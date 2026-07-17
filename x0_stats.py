#!/usr/bin/env python3
# x0_stats.py -- posterior stats of the exposure X0 = (1/2) ln(H0_local/H0)
# from the Model B production chain. Run from ~/class_holo_test:
#   python3 x0_stats.py
import glob, numpy as np

PREFIX = "chains_test/modelB_A2_v5_acc"
BURN = 0.3

files = sorted(f for f in glob.glob(PREFIX + ".*.txt")
               if f.split(".")[-2].isdigit())
assert files, f"no chain files found for {PREFIX}"

# column names from the cobaya header of the first file
with open(files[0]) as fh:
    header = fh.readline().lstrip("#").split()
iH0, iHloc, iw = header.index("H0"), header.index("H0_local"), 0

rows = []
for f in files:
    d = np.loadtxt(f)
    if d.ndim == 2 and len(d) > 10:
        rows.append(d[int(BURN * len(d)):])
d = np.vstack(rows)
w = d[:, iw]
H0, Hloc = d[:, iH0], d[:, iHloc]
X0 = 0.5 * np.log(Hloc / H0)

def wmean(x): return np.sum(w * x) / np.sum(w)
def wstd(x):  m = wmean(x); return np.sqrt(wmean((x - m) ** 2))
def wcorr(x, y):
    mx, my = wmean(x), wmean(y)
    return wmean((x - mx) * (y - my)) / (wstd(x) * wstd(y))

print(f"samples used        : {len(X0)}  (burn {BURN:.0%}, {len(files)} chains)")
print(f"<X0>                : {wmean(X0):.5f}")
print(f"sigma(X0)           : {wstd(X0):.5f}")
print(f"Corr(X0, H0_phys)   : {wcorr(X0, H0):+.3f}")
print(f"[check] e^(2<X0>)   : {np.exp(2 * wmean(X0)):.5f}   (expect ~1.0744)")
print(f"[check] hard gate   : 0.03450 ; realised interacting ~0.0359")
