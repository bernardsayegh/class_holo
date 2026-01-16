#!/usr/bin/env python3
from pathlib import Path

p = Path("source/perturbations.c")
lines = p.read_text().splitlines(True)

# locate the HOLOGRAPHIC INTERACTION banner line (you showed it's around 9238)
start = None
for i, ln in enumerate(lines):
    if "HOLOGRAPHIC INTERACTION (PERTURBATIONS)" in ln:
        start = i
        break

if start is None:
    raise SystemExit("[FATAL] Couldn't find 'HOLOGRAPHIC INTERACTION (PERTURBATIONS)' in source/perturbations.c")

# If banner line closes the comment on the same line (contains */), remove the closing so it opens a block
if "*/" in lines[start]:
    lines[start] = lines[start].replace("*/", "")  # keep opening /*, remove premature close

# Now find first "real code" line after the banner/doc lines:
# we treat these as "doc lines": blank, or starting with whitespace + '*', or starting with '/*'
end_insert = None
for j in range(start + 1, len(lines)):
    s = lines[j].lstrip()
    if s == "" or s == "\n":
        continue
    if s.startswith("*") or s.startswith("/*"):
        continue
    # found first line that isn't doc/comment-style -> insert closing */ before it
    end_insert = j
    break

if end_insert is None:
    raise SystemExit("[FATAL] Couldn't find where the doc block ends (unexpected file structure).")

# Only insert closing if there's no closing already between start and end_insert
chunk = "".join(lines[start:end_insert])
if "*/" not in chunk:
    lines.insert(end_insert, "*/\n")

p.write_text("".join(lines))
print(f"[OK] Converted holographic doc block into a proper /* ... */ comment (around line {start+1}).")
