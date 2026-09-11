#!/usr/bin/env python3
"""
Flag monoexonic transcripts that are on the wrong strand and flip them.

A transcript is considered misoriented when:
  - it is monoexonic, AND
  - its single exon has NO overlap with any miniprot alignment on the same strand, AND
  - its single exon overlaps at least one miniprot alignment on the OPPOSITE strand.

All supplied assembly files (StringTie GTF, Aletsch GTF, Trinity GFF3) are
corrected in place: the strand field of every line belonging to a misoriented
transcript is flipped.  Empty files (NO_FILE sentinels) are passed through
unchanged.

Usage:
    correct_misoriented.py \\
        --miniprot miniprot.gtf \\
        --input  stringtie.gtf aletsch.gtf trinity.gff \\
        --output corrected_stringtie.gtf corrected_aletsch.gtf corrected_trinity.gff
"""
import argparse
import os
import sys
from collections import defaultdict


def flip(strand):
    return '-' if strand == '+' else '+'


# ── attribute parsers ──────────────────────────────────────────────────────────

def tid_from_gtf(attrs):
    """Return transcript_id value from a GTF attributes string."""
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith('transcript_id'):
            tokens = part.split(None, 1)
            if len(tokens) == 2:
                return tokens[1].strip().strip('"')
    return None


def id_parent_from_gff3(attrs):
    """Return (ID, Parent) from a GFF3 attributes string."""
    id_ = parent = None
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith('ID='):
            id_ = part[3:]
        elif part.startswith('Parent='):
            parent = part[7:]
    return id_, parent


# ── helpers ────────────────────────────────────────────────────────────────────

def overlaps(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1


def any_overlap(chrom, strand, start, end, index):
    for s, e in index.get(chrom, {}).get(strand, []):
        if overlaps(start, end, s, e):
            return True
    return False


# ── parsing ────────────────────────────────────────────────────────────────────

def load_miniprot(path):
    """Return dict[chrom][strand] = [(start, end)] from miniprot GTF."""
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('transcript', 'mRNA'):
                continue
            idx[f[0]][f[6]].append((int(f[3]), int(f[4])))
    return idx


def collect_misoriented(input_files, miniprot_idx):
    """Return set of transcript IDs whose strand should be flipped."""
    misoriented = set()

    for path in input_files:
        if not os.path.getsize(path):
            continue

        exons     = defaultdict(list)   # tid -> [(start, end)]
        t_strand  = {}                  # tid -> strand
        t_chrom   = {}                  # tid -> chrom

        with open(path) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or f[2] != 'exon':
                    continue
                chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
                attrs = f[8]

                if 'transcript_id' in attrs:
                    tid = tid_from_gtf(attrs)
                else:
                    _, tid = id_parent_from_gff3(attrs)   # Parent= for exons

                if tid:
                    exons[tid].append((start, end))
                    t_strand[tid] = strand
                    t_chrom[tid]  = chrom

        for tid, ex in exons.items():
            if len(ex) != 1:
                continue
            chrom  = t_chrom[tid]
            strand = t_strand[tid]
            s, e   = ex[0]
            if (any_overlap(chrom, flip(strand), s, e, miniprot_idx)
                    and not any_overlap(chrom, strand, s, e, miniprot_idx)):
                misoriented.add(tid)

    return misoriented


# ── correction ─────────────────────────────────────────────────────────────────

# Feature types that carry a meaningful strand (gene lines are intentionally
# excluded; Mikado rebuilds gene coordinates from transcripts anyway)
_STRAND_FEATS = frozenset({'mRNA', 'transcript', 'exon', 'CDS',
                            'start_codon', 'stop_codon', 'three_prime_UTR',
                            'five_prime_UTR', 'UTR'})


def correct_file(in_path, out_path, misoriented):
    if not os.path.getsize(in_path):
        open(out_path, 'w').close()
        return

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

            # Determine the transcript ID for this line
            if 'transcript_id' in attrs:
                tid = tid_from_gtf(attrs)
            else:
                id_, parent = id_parent_from_gff3(attrs)
                # mRNA/transcript: use ID; children: use Parent
                tid = id_ if feat in ('mRNA', 'transcript') else parent

            if tid in misoriented and feat in _STRAND_FEATS:
                f[6] = flip(f[6])

            out.write('\t'.join(f) + '\n')


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--miniprot', required=True,
                    help='miniprot GTF (output of: miniprot --gtf)')
    ap.add_argument('--input',  nargs='+', required=True,
                    help='Assembly GFF/GTF files to correct')
    ap.add_argument('--output', nargs='+', required=True,
                    help='Output paths (one per --input file, same order)')
    args = ap.parse_args()

    if len(args.input) != len(args.output):
        sys.exit('--input and --output must have the same number of files')

    real_inputs = [p for p in args.input if os.path.getsize(p)]

    miniprot_idx = load_miniprot(args.miniprot)
    misoriented  = collect_misoriented(real_inputs, miniprot_idx)

    print(f'correct_misoriented: flipping {len(misoriented)} monoexonic transcripts',
          flush=True)

    for inp, out in zip(args.input, args.output):
        correct_file(inp, out, misoriented)


if __name__ == '__main__':
    main()
