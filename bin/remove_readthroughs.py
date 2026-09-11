#!/usr/bin/env python3
"""
Remove readthrough transcripts from a Mikado pick GFF3.

A readthrough is an mRNA whose exons overlap two or more distinct,
non-overlapping protein loci (as defined by merged miniprot alignments on
the same strand).  When ALL mRNAs of a gene are readthroughs, the gene
record is removed too.

Usage:
    remove_readthroughs.py \\
        --mikado  mikado_pick.loci.gff3 \\
        --miniprot miniprot.gtf \\
        --out     readthrough_filtered.gff3 \\
        --report  readthroughs.txt
"""
import argparse
import sys
from collections import defaultdict


# ── interval utilities ─────────────────────────────────────────────────────────

def merge_intervals(ivs):
    """Return sorted, non-overlapping merged list of (start, end) intervals."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


# ── parsing ────────────────────────────────────────────────────────────────────

def get_attr(attrs, key):
    """Return value for key= in a semicolon-separated GFF3 attribute string."""
    for part in attrs.split(';'):
        part = part.strip()
        if part.startswith(key + '='):
            return part[len(key) + 1:]
    return None


def load_protein_loci(path):
    """
    Parse miniprot GTF; return merged protein loci.
    Returns dict[chrom][strand] = [(start, end), ...] sorted and merged.
    """
    raw = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('transcript', 'mRNA'):
                continue
            raw[f[0]][f[6]].append((int(f[3]), int(f[4])))
    loci = {}
    for chrom, strands in raw.items():
        loci[chrom] = {}
        for strand, ivs in strands.items():
            loci[chrom][strand] = merge_intervals(ivs)
    return loci


def loci_hit_by_exons(exons, loci_list):
    """Return set of locus indices (into loci_list) overlapped by any exon."""
    hit = set()
    for es, ee in exons:
        for i, (ls, le) in enumerate(loci_list):
            if es <= le and ls <= ee:
                hit.add(i)
    return hit


def find_readthroughs(mikado_path, protein_loci):
    """
    Return (readthrough_mRNA_ids, gene_mrna_map).
    gene_mrna_map: gene_id -> list of mRNA IDs.
    """
    # First pass: collect gene→mRNA mapping and mRNA exons
    gene_mrnas = defaultdict(list)   # gene_id -> [mrna_id, ...]
    mrna_exons = defaultdict(list)   # mrna_id -> [(start, end), ...]
    mrna_meta  = {}                   # mrna_id -> (chrom, strand)

    with open(mikado_path) as fh:
        cur_mrna = None
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            feat, attrs = f[2], f[8]
            if feat == 'mRNA':
                mid = get_attr(attrs, 'ID')
                pid = get_attr(attrs, 'Parent')
                cur_mrna = mid
                if mid and pid:
                    gene_mrnas[pid].append(mid)
                    mrna_meta[mid] = (f[0], f[6])
            elif feat == 'exon':
                parent = get_attr(attrs, 'Parent')
                if parent:
                    mrna_exons[parent].append((int(f[3]), int(f[4])))

    readthrough_ids = set()
    for mid, exons in mrna_exons.items():
        if mid not in mrna_meta:
            continue
        chrom, strand = mrna_meta[mid]
        if chrom not in protein_loci or strand not in protein_loci[chrom]:
            continue
        hit = loci_hit_by_exons(exons, protein_loci[chrom][strand])
        if len(hit) >= 2:
            readthrough_ids.add(mid)

    return readthrough_ids, dict(gene_mrnas)


# ── GFF3 filtering ─────────────────────────────────────────────────────────────

def write_filtered(mikado_path, out_path, readthrough_ids, gene_mrnas):
    """Stream the GFF3, dropping readthrough mRNAs and genes that lose all mRNAs."""
    genes_to_drop = {
        gid for gid, mrnas in gene_mrnas.items()
        if mrnas and all(m in readthrough_ids for m in mrnas)
    }

    skip_mrna = None
    with open(mikado_path) as fh, open(out_path, 'w') as out:
        for line in fh:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                out.write(line)
                continue
            feat, attrs = f[2], f[8]

            if feat in ('gene', 'ncRNA_gene'):
                skip_mrna = None
                gid = get_attr(attrs, 'ID')
                if gid in genes_to_drop:
                    continue
                out.write(line)

            elif feat == 'mRNA':
                mid = get_attr(attrs, 'ID')
                if mid in readthrough_ids:
                    skip_mrna = mid
                else:
                    skip_mrna = None
                    out.write(line)

            else:
                parent = get_attr(attrs, 'Parent')
                if parent == skip_mrna:
                    continue
                out.write(line)


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mikado',   required=True, help='Mikado pick GFF3')
    ap.add_argument('--miniprot', required=True, help='miniprot GTF')
    ap.add_argument('--out',      required=True, help='Output filtered GFF3')
    ap.add_argument('--report',   default=None,  help='Optional: write readthrough IDs here')
    args = ap.parse_args()

    protein_loci = load_protein_loci(args.miniprot)
    readthrough_ids, gene_mrnas = find_readthroughs(args.mikado, protein_loci)

    print(f'remove_readthroughs: found {len(readthrough_ids)} readthrough transcript(s)',
          flush=True)

    if args.report:
        with open(args.report, 'w') as rp:
            for mid in sorted(readthrough_ids):
                rp.write(mid + '\n')

    write_filtered(args.mikado, args.out, readthrough_ids, gene_mrnas)


if __name__ == '__main__':
    main()
