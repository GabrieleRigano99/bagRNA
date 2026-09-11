#!/usr/bin/env python3
"""Convert compleasm protein full_table.tsv + scores_cutoff → Mikado external busco_score TSV.

Protein IDs in TD2 longest_orfs.pep are TranscriptID.pN; extract transcript by stripping .pN.
For each transcript, keep the highest score/cutoff ratio across all BUSCO hits (capped at 1.0).
Only outputs transcripts with at least one passing (Single/Duplicated) BUSCO hit.
"""
import re
import sys
from collections import defaultdict


def parse_scores_cutoff(path):
    """Return dict keyed by both full ID (e.g. '164661at147550') and plain ID ('164661')."""
    cutoffs = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                val = float(parts[1])
                cutoffs[parts[0]] = val
                cutoffs[parts[0].split("at")[0]] = val  # plain ID fallback
    return cutoffs


def protein_to_transcript(protein_id):
    return re.sub(r'\.p\d+$', '', protein_id)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} full_table.tsv scores_cutoff")

    full_table_path = sys.argv[1]
    scores_cutoff_path = sys.argv[2]

    cutoffs = parse_scores_cutoff(scores_cutoff_path)

    best = defaultdict(float)  # transcript → best busco_score ratio

    with open(full_table_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            busco_full = parts[0]
            # Strip "atNNNNNN" lineage suffix if present: e.g. "157615at147550" → "157615"
            busco_id = busco_full.split("at")[0]
            status = parts[1]
            protein_id = parts[2]
            try:
                score = float(parts[3])
            except ValueError:
                continue
            if status not in ("Single", "Duplicated", "Fragmented"):
                continue
            cutoff = cutoffs.get(busco_id) or cutoffs.get(busco_full)
            if cutoff is None or cutoff == 0:
                continue
            transcript = protein_to_transcript(protein_id)
            ratio = min(1.0, score / cutoff)
            if ratio > best[transcript]:
                best[transcript] = ratio

    print("tid\tbusco_score")
    for tid, score in sorted(best.items()):
        print(f"{tid}\t{score:.6f}")


if __name__ == "__main__":
    main()
