#!/usr/bin/env python3
"""
Fix TD2 v1.0.7/v1.1.0 BED output: minus-strand ORFs produce an empty thickEnd
because the '-' inside '(-)' breaks the naive split used in Predict.py.
Re-extracts the end coordinate from the name field via regex.
Usage: fix_td2_bed.py <bed_file>  (edits in-place)
"""
import re
import os
import sys

bed = sys.argv[1]
fixed = bed + ".fix"

with open(bed) as fin, open(fixed, "w") as fout:
    for line in fin:
        if line.startswith("track") or "\t" not in line:
            fout.write(line)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 8 and f[7] == "":
            m = re.search(r":(\d+)-(\d+)\([+-]\)$", f[3])
            if m:
                f[7] = m.group(2)
        fout.write("\t".join(f) + "\n")

os.replace(fixed, bed)
