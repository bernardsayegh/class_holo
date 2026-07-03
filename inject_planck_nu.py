#!/usr/bin/env python3
"""Inject the Planck-standard neutrino convention (v5 publication grid)
into every configs/*.yaml. Dry-run by default; use --apply to write.

    python3 inject_planck_nu.py           # show what would change
    python3 inject_planck_nu.py --apply   # write changes
"""
import glob, re, sys

APPLY = "--apply" in sys.argv
BLOCK = ["N_ncdm: 1", "m_ncdm: 0.06", "N_ur: 2.0328"]

changed, skipped, failed = [], [], []
for fn in sorted(glob.glob("configs/*.yaml")):
    s = open(fn).read()
    if "m_ncdm" in s:
        skipped.append(fn); continue
    m = list(re.finditer(r"^(\s*)extra_args:\s*$", s, re.M))
    if len(m) != 1:
        failed.append((fn, f"{len(m)} extra_args blocks")); continue
    indent = m[0].group(1) + "  "
    insert = "".join(f"{indent}{line}\n" for line in BLOCK)
    pos = m[0].end() + 1
    new = s[:pos] + insert + s[pos:]
    if APPLY:
        open(fn, "w").write(new)
    changed.append(fn)

print(f"{'APPLIED' if APPLY else 'DRY RUN'}:")
for f in changed: print(f"  + nu block -> {f}")
for f in skipped: print(f"  = already has m_ncdm: {f}")
for f, why in failed: print(f"  ! MANUAL FIX NEEDED ({why}): {f}")
if failed: sys.exit(1)
print(f"\n{len(changed)} injected, {len(skipped)} skipped, 0 failures.")
if not APPLY: print("Re-run with --apply to write.")
