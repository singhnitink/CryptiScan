#!/usr/bin/env python3
"""
check_stage10.py -- build-time regression guard for guided pocket detection.

Called from the container's %test block. Every assertion here corresponds to a
bug that made Stage 10 return "no pocket" for every frame while looking like a
genuine negative result. Each one cost a full 2 x 2000-frame run to find, so
the build fails rather than let one return silently.

Run standalone to check a working tree:
    python3 check_stage10.py [path/to/cryptic_pocket_analysis.py]
"""

import re
import sys
from pathlib import Path

TARGET = "/opt/pipeline/cryptic_pocket_analysis.py"

CHECKS = [
    (r"\{r\}:-:\{frame_chain\}",
     "-P residue spec must be resnum:-:chain.",
     "fpocket splits each residue on ':' with strtok_r, which collapses an "
     "empty field, and the switch consuming the pieces has no break "
     "statements. With '127::A' the insertion code ends up holding the chain "
     "letter and the match fails for every residue, so nothing is ever found."),

    (r'"-u",\s*str\(min_explicit\)',
     "-u must be passed to fpocket.",
     "An alpha sphere is kept only if at least min_n_explicit_pocket_atoms of "
     "its four contact atoms belong to the listed residues. fpocket's default "
     "is 4, i.e. all four, which for a site of five residues retains almost "
     "nothing: every frame returns 'No Pockets Found while refining'."),

    (r"pocket\{pnum\}_atm\.pdb",
     "Pocket atom files are numbered from 1, matching info.txt.",
     "'Pocket 1 :' in info.txt corresponds to pockets/pocket1_atm.pdb. "
     "Subtracting one looks for pocket0_atm.pdb, which never exists, so every "
     "pocket is skipped and every frame is reported closed."),

    (r"safe_label",
     "The per-ensemble temp directory name must be sanitised.",
     "fpocket builds its output paths with C string operations and cannot "
     "write into a directory whose name contains a space. The label is "
     "'wild type', so an unsanitised prefix makes EVERY wild-type frame fail "
     "while the mutant succeeds, which reads as a one-sided biological "
     "result rather than a path bug."),

    (r"n_pockets_unguided",
     "The specificity control needs its own unguided pass.",
     "A guided run reports only the cryptic pocket, so summing 'every other "
     "pocket' over it is always zero. Without a separate unguided pass the "
     "control that separates a site-specific opening from a general "
     "loosening of the surface silently reports no background at all."),

    (r"rmsd_ownref",
     "Backbone RMSD must use one common reference for both ensembles.",
     "Referring each ensemble to its own first frame makes the two values "
     "incomparable: it produced a large apparent effect in a control that "
     "should be flat, purely because one frame 0 sat further from its own "
     "centroid."),
]


def main():
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(TARGET)
    if not target.is_file():
        print(f"  FAIL    {target} not found")
        return 1

    src = target.read_text()
    failed = []
    for pattern, headline, why in CHECKS:
        if not re.search(pattern, src):
            failed.append((headline, why))

    for headline, why in failed:
        print(f"  FAIL    {headline}")
        for line in (why[i:i + 68] for i in range(0, len(why), 68)):
            print(f"          {line}")

    if failed:
        print(f"\n  {len(failed)} of {len(CHECKS)} guards failed in "
              f"{target.name}.")
        return 1

    print(f"  all {len(CHECKS)} Stage 10 guards present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
