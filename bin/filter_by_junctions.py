#!/usr/bin/env python3
"""
Filter assembled RNAseq transcripts to keep only those whose splice junctions
are all validated by Portcullis.

Rules:
  - Monoexonic transcripts: always kept (no junctions to validate).
  - Multi-exonic transcripts: kept only when EVERY consecutive exon-pair
    produces a junction key (chrom, exon_n_end, exon_n+1_start-1) that is
    present in the Portcullis BED12 pass set.
  - Empty / NO_FILE input files are passed through unchanged.

Portcullis BED12 column convention (0-based field index):
  0: chrom   5: strand   6: thickStart   7: thickEnd
  thickStart = first intron base (= GFF exon1_end, 1-based)
  thickEnd   = last intron base + 1 exclusive (= GFF exon2_start - 1, 1-based)

Usage:
    filter_by_junctions.py \\
        --junctions portcullis.pass.junctions.bed \\
        --input  stringtie.gtf aletsch.gtf trinity.gff \\
        --output validated_stringtie.gtf validated_aletsch.gtf validated_trinity.gff
"""
import argparse
import os
import sys
from collections import defaultdict


# ── Portcullis BED12 parsing ───────────────────────────────────────────────────

def load_junctions(bed_path):
    """
    Return a frozenset of (chrom, thickStart, thickEnd) tuples from a
    Portcullis BED12 pass file.  Strand is intentionally ignored so the
    set can be used against strand-corrected assembled transcripts without
    requiring exact strand agreement.
    """
    junctions = set()
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('track') or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            try:
                junctions.add((f[0], int(f[6]), int(f[7])))
            except ValueError:
                continue
    return frozenset(junctions)


# ── attribute helpers ─────────────────────────────────────────────────────────

def tid_from_gtf(attrs):
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith('transcript_id'):
            tokens = part.split(None, 1)
            return tokens[1].strip().strip('"') if len(tokens) == 2 else None
    return None


def tid_from_gff3_exon(attrs):
    """Return Parent= value from a GFF3 exon attribute string."""
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith('Parent='):
            return part[7:]
    return None


def id_from_gff3(attrs):
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith('ID='):
            return part[3:]
    return None


# ── per-file filtering ────────────────────────────────────────────────────────

def collect_invalid_transcripts(path, junctions):
    """
    Return set of transcript IDs that should be removed because they have
    at least one multi-exon junction absent from the Portcullis set.
    """
    # Collect exons per transcript
    exons    = defaultdict(list)   # tid -> [(chrom, start, end)]
    is_gtf   = False

    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                if 'transcript_id' in (f[8] if len(f) >= 9 else ''):
                    is_gtf = True
                continue

            chrom, start, end = f[0], int(f[3]), int(f[4])
            attrs = f[8]

            if 'transcript_id' in attrs:
                is_gtf = True
                tid = tid_from_gtf(attrs)
            else:
                tid = tid_from_gff3_exon(attrs)

            if tid:
                exons[tid].append((chrom, start, end))

    invalid = set()
    for tid, ex_list in exons.items():
        if len(ex_list) <= 1:
            continue   # monoexonic → always keep
        # Sort exons by start position to get correct consecutive pairs
        ex_sorted = sorted(ex_list, key=lambda x: x[1])
        for i in range(len(ex_sorted) - 1):
            chrom  = ex_sorted[i][0]
            e1_end = ex_sorted[i][2]       # exon_n end (1-based)
            e2_beg = ex_sorted[i + 1][1]   # exon_n+1 start (1-based)
            key    = (chrom, e1_end, e2_beg - 1)
            if key not in junctions:
                invalid.add(tid)
                break

    return invalid


def write_filtered_gtf(in_path, out_path, invalid):
    """Stream GTF removing transcript and exon lines for invalid transcript IDs."""
    skip_tid = None
    with open(in_path) as fh, open(out_path, 'w') as out:
        for line in fh:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                out.write(line)
                continue
            feat, attrs = f[2], f[8]

            if feat == 'gene':
                # Gene records span multiple transcripts; don't skip by gene
                skip_tid = None
                out.write(line)
                continue

            tid = tid_from_gtf(attrs)
            if tid is None:
                out.write(line)
                continue

            if tid in invalid:
                skip_tid = tid
                continue   # skip this line (transcript or exon)

            if tid != skip_tid:
                skip_tid = None
                out.write(line)

    return len(invalid)


def write_filtered_gff3(in_path, out_path, invalid):
    """Stream GFF3 removing mRNA and child lines for invalid transcript IDs."""
    skip_mid = None
    with open(in_path) as fh, open(out_path, 'w') as out:
        for line in fh:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                out.write(line)
                continue
            feat, attrs = f[2], f[8]

            if feat in ('gene',):
                skip_mid = None
                out.write(line)
                continue

            if feat in ('mRNA', 'transcript'):
                mid = id_from_gff3(attrs)
                if mid in invalid:
                    skip_mid = mid
                    continue
                skip_mid = None
                out.write(line)
                continue

            # Child feature: skip if its Parent is the skipped mRNA
            parent = tid_from_gff3_exon(attrs)  # reads Parent=
            if parent == skip_mid:
                continue
            out.write(line)


def filter_file(in_path, out_path, junctions):
    if not os.path.getsize(in_path):
        open(out_path, 'w').close()
        return

    invalid = collect_invalid_transcripts(in_path, junctions)
    print(f'  {os.path.basename(in_path)}: removing {len(invalid)} '
          f'transcript(s) with unvalidated junctions', flush=True)

    # Detect format by scanning a few non-comment lines
    is_gtf = False
    with open(in_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.split('\t')
            if len(f) >= 9 and 'transcript_id' in f[8]:
                is_gtf = True
            break

    if is_gtf:
        write_filtered_gtf(in_path, out_path, invalid)
    else:
        write_filtered_gff3(in_path, out_path, invalid)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--junctions', required=True,
                    help='Portcullis filtered pass BED12 file')
    ap.add_argument('--input',  nargs='+', required=True)
    ap.add_argument('--output', nargs='+', required=True)
    args = ap.parse_args()

    if len(args.input) != len(args.output):
        sys.exit('--input and --output must have the same number of files')

    junctions = load_junctions(args.junctions)
    print(f'filter_by_junctions: loaded {len(junctions)} validated junctions',
          flush=True)

    for inp, out in zip(args.input, args.output):
        filter_file(inp, out, junctions)


if __name__ == '__main__':
    main()
