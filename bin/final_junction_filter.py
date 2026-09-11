#!/usr/bin/env python3
"""Final junction filter: remove multi-exonic gene models with zero portcullis support.

Rule per gene:
  - A model is "supported" if: monoexonic, ProteinEvidence, or has at least one
    portcullis-confirmed splice junction.
  - If the gene has at least one supported model:
      Drop any unsupported models (zero confirmed junctions).
  - If ALL models are unsupported:
      Keep them all as-is — a multi-exonic model with unconfirmed junctions
      produces a better protein than a collapsed monoexonic one.

Usage:
    final_junction_filter.py input.gff3 portcullis.bed[,...] output.gff3
"""
import sys
from collections import defaultdict


TRANSCRIPT_TYPES = frozenset(['mRNA', 'ncRNA', 'lncRNA', 'transcript'])


def attr(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


def load_portcullis(paths):
    junctions = set()
    for path in paths:
        with open(path) as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 8:
                    continue
                try:
                    junctions.add((cols[0], int(cols[6]), int(cols[7])))
                except (ValueError, IndexError):
                    continue
    return frozenset(junctions)


def junctions_of(exons, chrom):
    s = sorted(exons)
    return [(chrom, s[i][1], s[i + 1][0] - 1) for i in range(len(s) - 1)]



def main():
    args = sys.argv[1:]
    if len(args) < 3:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 portcullis.bed[,...] output.gff3")

    in_path   = args[0]
    out_path  = args[-1]
    bed_paths = args[1:-1]

    portcullis = load_portcullis(bed_paths)
    sys.stderr.write(f"Loaded {len(portcullis)} portcullis junctions\n")

    with open(in_path) as f:
        raw = f.readlines()

    # ── First pass: collect structure ─────────────────────────────────────────
    gene_tids = defaultdict(list)   # gid → [tid, ...]  in order
    mrnas     = {}                  # tid → {gene, chrom, strand, start, end, exons, cds, alias}
    cur_tid   = None
    known_genes = set()

    for line in raw:
        if line.startswith('#') or not line.strip():
            continue
        c = line.split('\t')
        if len(c) < 9:
            continue
        a    = attr(c[8])
        feat = c[2]

        if feat == 'gene':
            known_genes.add(a.get('ID', ''))
            cur_tid = None

        elif feat in TRANSCRIPT_TYPES:
            parent = a.get('Parent', '')
            if parent not in known_genes:
                cur_tid = None
                continue
            cur_tid = a.get('ID', '')
            mrnas[cur_tid] = {
                'gene':   parent,
                'chrom':  c[0], 'strand': c[6],
                'start':  int(c[3]), 'end': int(c[4]),
                'exons':  [], 'cds': [],
                'alias':  a.get('alias', ''),
            }
            gene_tids[parent].append(cur_tid)

        elif feat == 'exon':
            parent = a.get('Parent', cur_tid)
            if parent in mrnas:
                mrnas[parent]['exons'].append((int(c[3]), int(c[4])))

        elif feat == 'CDS':
            parent = a.get('Parent', cur_tid)
            if parent in mrnas:
                mrnas[parent]['cds'].append((int(c[3]), int(c[4])))

    # ── Classify mRNAs and decide fate ───────────────────────────────────────
    def is_supported(tid):
        m = mrnas[tid]
        alias = m['alias'].lower()
        if 'proteinevidence' in alias:
            return True                           # exempt
        if len(m['exons']) <= 1:
            return True                           # monoexonic — no junctions to verify
        js = junctions_of(m['exons'], m['chrom'])
        return any(j in portcullis for j in js)  # at least one junction confirmed

    drop_tid = set()
    n_dropped = 0

    for gid, tids in gene_tids.items():
        supported   = [t for t in tids if is_supported(t)]
        unsupported = [t for t in tids if not is_supported(t)]

        if not unsupported:
            continue

        if supported:
            # Good models available — drop models with zero confirmed junctions
            for t in unsupported:
                drop_tid.add(t)
                n_dropped += 1
        # else: all models unsupported — keep them all as-is for BUSCO coverage

    sys.stderr.write(f"Dropped (zero portcullis support, alternatives available): {n_dropped}\n")

    # ── Second pass: write corrected GFF ─────────────────────────────────────
    skip_parent = None

    with open(out_path, 'w') as out:
        for line in raw:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                out.write(line)
                continue
            feat = c[2]
            a    = attr(c[8])

            if feat == 'gene':
                skip_parent = None
                out.write(line)

            elif feat in TRANSCRIPT_TYPES:
                tid = a.get('ID', '')
                if tid in drop_tid:
                    skip_parent = tid
                    continue
                skip_parent = None
                out.write(line)

            else:
                if a.get('Parent', '') == skip_parent:
                    continue
                out.write(line)


if __name__ == '__main__':
    main()
