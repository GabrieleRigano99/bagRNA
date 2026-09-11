#!/usr/bin/env python3
"""Add gene loci missing from the current annotation from two sources:
  1. TD2 genome GFF3 (TransDecoder predictions on mikado_prepared transcripts):
     Add any gene that does not overlap an existing gene on the same strand.
  2. Compleasm genome GFF (miniprot predictions from BUSCO proteins):
     Group overlapping mRNAs into loci, pick best (highest score),
     add loci that still don't overlap any existing or newly-added gene.

Usage:
    add_missing_loci.py current.gff3 td2_genome.gff3 busco.gff output.gff3
"""
import sys
import bisect
from collections import defaultdict


def attr(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


# ── Interval index ─────────────────────────────────────────────────────────────
class IntervalIndex:
    """Per (chrom, strand): sorted list of (start, end) for fast overlap queries."""

    def __init__(self):
        self._starts = defaultdict(list)   # key → sorted starts
        self._ends   = defaultdict(list)   # key → ends parallel to _starts

    def add(self, chrom, strand, start, end):
        key = (chrom, strand)
        i = bisect.bisect_left(self._starts[key], start)
        self._starts[key].insert(i, start)
        self._ends[key].insert(i, end)

    def overlaps(self, chrom, strand, start, end):
        key = (chrom, strand)
        starts = self._starts[key]
        ends   = self._ends[key]
        if not starts:
            return False
        # All intervals with start <= end
        hi = bisect.bisect_right(starts, end)
        for i in range(hi):
            if ends[i] >= start:
                return True
        return False


# ── Parse current annotation ───────────────────────────────────────────────────
def load_current(path):
    idx = IntervalIndex()
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.split('\t')
            if len(c) < 9 or c[2] != 'gene':
                continue
            a = attr(c[8])
            idx.add(c[0], c[6], int(c[3]), int(c[4]))
    return idx


# ── TD2 genome GFF3 ────────────────────────────────────────────────────────────
def add_td2_loci(td2_path, idx, out):
    """Stream TD2 genome GFF3; emit gene blocks that don't overlap any indexed gene."""
    n_added = n_skipped = 0
    cur_gene = None
    cur_chrom = cur_strand = None
    cur_start = cur_end = 0
    buf = []

    def flush():
        nonlocal n_added, n_skipped
        if cur_gene is None:
            return
        if idx.overlaps(cur_chrom, cur_strand, cur_start, cur_end):
            n_skipped += 1
        else:
            for l in buf:
                out.write(l)
            idx.add(cur_chrom, cur_strand, cur_start, cur_end)
            n_added += 1

    with open(td2_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.split('\t')
            if len(c) < 9:
                continue
            if c[2] == 'gene':
                flush()
                buf = [line]
                a = attr(c[8])
                cur_gene  = a.get('ID', '')
                cur_chrom = c[0]
                cur_strand = c[6]
                cur_start  = int(c[3])
                cur_end    = int(c[4])
            else:
                buf.append(line)

    flush()
    sys.stderr.write(f"TD2 loci: added={n_added}, skipped (overlap)={n_skipped}\n")


# ── Compleasm busco.gff ────────────────────────────────────────────────────────
def add_busco_loci(busco_path, idx, out):
    """
    Group overlapping mRNAs on the same (chrom, strand) into loci.
    Pick the best mRNA (highest score) per locus.
    Add loci that don't overlap any indexed gene.
    Synthesise gene + exon features from mRNA and CDS records.
    """
    # Collect all mRNAs
    mrnas = []   # (chrom, strand, start, end, score, mrna_id, [cds_lines])
    cur_mrna = None
    cds_buf  = []

    with open(busco_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            feat = c[2]
            if feat == 'mRNA':
                if cur_mrna is not None:
                    mrnas.append(cur_mrna + (cds_buf,))
                a = attr(c[8])
                try:
                    score = float(c[5])
                except ValueError:
                    score = 0.0
                cur_mrna = (c[0], c[6], int(c[3]), int(c[4]), score, a.get('ID', ''))
                cds_buf  = [line]
            elif feat in ('CDS', 'stop_codon') and cur_mrna is not None:
                cds_buf.append(line)

    if cur_mrna is not None:
        mrnas.append(cur_mrna + (cds_buf,))

    # Sort by chrom, strand, start
    mrnas.sort(key=lambda x: (x[0], x[1], x[2]))

    # Cluster into loci
    loci = []   # each: [(chrom, strand, start, end, score, mrna_id, [lines])]
    for m in mrnas:
        chrom, strand, start, end = m[0], m[1], m[2], m[3]
        if loci and loci[-1][0][0] == chrom and loci[-1][0][1] == strand and start <= loci[-1][0][3]:
            # overlaps last locus — extend and add
            prev = loci[-1]
            # update locus end if needed
            loci[-1] = prev
            prev.append(m)
            # extend locus span
            loci[-1][0] = (chrom, strand, prev[0][2], max(prev[0][3], end)) + prev[0][4:]
        else:
            loci.append([m])

    # Simpler clustering: separate pass
    loci = []
    cur_locus = None
    cur_locus_end = 0

    for m in mrnas:
        chrom, strand, start, end = m[0], m[1], m[2], m[3]
        if (cur_locus is None
                or cur_locus[0][0] != chrom
                or cur_locus[0][1] != strand
                or start > cur_locus_end):
            cur_locus = [m]
            cur_locus_end = end
            loci.append(cur_locus)
        else:
            cur_locus.append(m)
            if end > cur_locus_end:
                cur_locus_end = end

    n_added = n_skipped = 0
    gene_counter = [0]

    for locus in loci:
        # locus span
        chrom   = locus[0][0]
        strand  = locus[0][1]
        l_start = min(m[2] for m in locus)
        l_end   = max(m[3] for m in locus)

        if idx.overlaps(chrom, strand, l_start, l_end):
            n_skipped += 1
            continue

        # Pick best mRNA by score
        best = max(locus, key=lambda m: m[4])
        chrom, strand, start, end, score, mrna_id, lines = best

        gene_counter[0] += 1
        gid = f"busco_gene_{gene_counter[0]:06d}"
        tid = f"busco_mrna_{gene_counter[0]:06d}"

        # Synthesise gene feature
        out.write(f"{chrom}\tcompleasm\tgene\t{start}\t{end}\t{score}\t{strand}\t.\t"
                  f"ID={gid};Name={mrna_id}\n")

        # mRNA feature
        out.write(f"{chrom}\tcompleasm\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t"
                  f"ID={tid};Parent={gid};Name={mrna_id}\n")

        # Collect CDS segments to synthesise exons
        cds_segments = []
        for l in lines:
            c = l.rstrip('\n').split('\t')
            if len(c) >= 9 and c[2] == 'CDS':
                cds_segments.append((int(c[3]), int(c[4]), c[7], l))

        # Exon = one per CDS segment (CDS-only gene, no UTR)
        for i, (cs, ce, phase, orig_line) in enumerate(sorted(cds_segments)):
            out.write(f"{chrom}\tcompleasm\texon\t{cs}\t{ce}\t.\t{strand}\t.\t"
                      f"ID={tid}.exon{i+1};Parent={tid}\n")

        # CDS features (re-emit with updated Parent)
        for i, (cs, ce, phase, orig_line) in enumerate(sorted(cds_segments)):
            out.write(f"{chrom}\tcompleasm\tCDS\t{cs}\t{ce}\t.\t{strand}\t{phase}\t"
                      f"ID=cds.{tid}.{i+1};Parent={tid}\n")

        idx.add(chrom, strand, start, end)
        n_added += 1

    sys.stderr.write(f"Busco loci: added={n_added}, skipped (overlap)={n_skipped}\n")


def main():
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} current.gff3 td2_genome.gff3 busco.gff output.gff3")

    current_path = sys.argv[1]
    td2_path     = sys.argv[2]
    busco_path   = sys.argv[3]
    out_path     = sys.argv[4]

    sys.stderr.write("Loading current annotation intervals...\n")
    idx = load_current(current_path)

    with open(out_path, 'w') as out:
        # Copy current annotation verbatim
        with open(current_path) as f:
            for line in f:
                out.write(line)

        out.write("##\n## --- TD2 non-overlapping loci ---\n##\n")
        sys.stderr.write("Adding TD2 loci...\n")
        add_td2_loci(td2_path, idx, out)

        out.write("##\n## --- Compleasm/miniprot BUSCO loci ---\n##\n")
        sys.stderr.write("Adding compleasm/BUSCO loci...\n")
        add_busco_loci(busco_path, idx, out)

    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
