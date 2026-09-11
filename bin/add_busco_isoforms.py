#!/usr/bin/env python3
"""For each BUSCO missing from the current protein annotation:
  1. Find the best miniprot model for that BUSCO in compleasm busco.gff.
  2. Find the gene in the current annotation that overlaps it (same chrom, any strand).
  3. Add the miniprot mRNA as an extra isoform; expand gene boundaries if needed.
  4. If no gene overlaps, add as a new standalone gene.

Usage:
    add_busco_isoforms.py current.gff3 busco.gff missing_ids.txt output.gff3
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


def attr_str(d):
    return ';'.join(f"{k}={v}" for k, v in d.items())


# ── Load missing BUSCO IDs ────────────────────────────────────────────────────
def load_missing(path):
    ids = set()
    with open(path) as f:
        for line in f:
            ids.add(line.strip())
    return ids


# ── Load busco.gff: collect best model per BUSCO ID ──────────────────────────
def load_busco_models(path, missing_ids):
    """Return best miniprot mRNA per BUSCO ID (highest score), only for missing IDs.
    Result: busco_id → (chrom, strand, start, end, score, mrna_id, [cds_lines])
    """
    # First pass: collect all mRNAs for missing BUSCOs
    # BUSCO ID in Target field: "287901at147550_..." → id = "287901"
    candidates = defaultdict(list)   # busco_id → [(chrom, strand, start, end, score, mrna_id, cds_buf)]
    cur_mrna = None
    cur_busco = None
    cds_buf = []

    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            feat = c[2]
            if feat == 'mRNA':
                if cur_mrna is not None and cur_busco in missing_ids:
                    candidates[cur_busco].append(
                        (cur_mrna[0], cur_mrna[1], cur_mrna[2], cur_mrna[3],
                         cur_mrna[4], cur_mrna[5], cds_buf)
                    )
                a = attr(c[8])
                target = a.get('Target', '')
                cur_busco = target.split('at')[0] if 'at' in target else ''
                try:
                    score = float(c[5])
                except ValueError:
                    score = 0.0
                cur_mrna = (c[0], c[6], int(c[3]), int(c[4]), score, a.get('ID', ''))
                cds_buf = []
            elif feat == 'CDS' and cur_mrna is not None:
                cds_buf.append(c)

    if cur_mrna is not None and cur_busco in missing_ids:
        candidates[cur_busco].append(
            (cur_mrna[0], cur_mrna[1], cur_mrna[2], cur_mrna[3],
             cur_mrna[4], cur_mrna[5], cds_buf)
        )

    # Return best model per BUSCO (highest score)
    best = {}
    for bid, mods in candidates.items():
        best[bid] = [max(mods, key=lambda m: m[4])]

    sys.stderr.write(f"Found miniprot models for {len(best)}/{len(missing_ids)} target BUSCOs\n")
    return best


# ── Load current annotation ───────────────────────────────────────────────────
def load_annotation(path):
    """Returns:
    - raw: all lines
    - gene_meta: gid → {chrom, strand, start, end, line_idx}
    - by_loc: (chrom, strand) → sorted [(start, end, gid)]
    - gene_order: [gid, ...] in file order
    - gene_cds_sigs: gid → set of CDS-signature tuples, one per existing mRNA
      (used to skip attaching a busco isoform that would duplicate an
      existing isoform's CDS — table2asn emits one CDS feature per mRNA,
      and antiSMASH rejects multiple CDS features at the same location)
    """
    raw = []
    with open(path) as f:
        raw = f.readlines()

    gene_meta = {}
    gene_order = []
    by_loc = defaultdict(list)
    gene_cds_sigs = defaultdict(set)
    cur_gid = None
    cur_tid = None
    cur_cds = []

    def flush_mrna():
        if cur_gid is not None and cur_cds:
            gene_cds_sigs[cur_gid].add(tuple(sorted(cur_cds)))

    for i, line in enumerate(raw):
        if line.startswith('#') or not line.strip():
            continue
        c = line.split('\t')
        if len(c) < 9:
            continue
        feat = c[2]
        a = attr(c[8])

        if feat == 'gene':
            flush_mrna()
            gid = a.get('ID', '')
            start, end = int(c[3]), int(c[4])
            gene_meta[gid] = {
                'chrom': c[0], 'strand': c[6],
                'start': start, 'end': end,
                'line_idx': i,
            }
            gene_order.append(gid)
            by_loc[(c[0], c[6])].append((start, end, gid))
            cur_gid = None
            cur_tid = None
            cur_cds = []
        elif feat in ('mRNA', 'ncRNA', 'lncRNA', 'transcript'):
            flush_mrna()
            cur_tid = a.get('ID', '')
            cur_gid = a.get('Parent', '')
            cur_cds = []
        elif feat == 'CDS':
            parent = a.get('Parent', cur_tid)
            if parent == cur_tid:
                cur_cds.append((int(c[3]), int(c[4])))
    flush_mrna()

    # Sort each locus list by start
    for key in by_loc:
        by_loc[key].sort()

    return raw, gene_meta, dict(by_loc), gene_order, dict(gene_cds_sigs)


def _near_dup_cds(a, b, tol=5):
    """True if two sorted CDS coordinate lists are the same isoform, differing
    only by a small (<=tol bp) offset at an outer boundary — e.g. a partial
    start/stop codon in a miniprot model. table2asn's frame-correction
    (-c ewf) can snap such near-identical CDS to the exact same final span,
    which antiSMASH then rejects as a duplicate CDS location.
    All internal splice junctions must match exactly; only the two outer
    boundaries may differ, and only within tol.
    """
    if len(a) != len(b):
        return False
    if len(a) == 1:
        (s1, e1), (s2, e2) = a[0], b[0]
        return abs(s1 - s2) <= tol and abs(e1 - e2) <= tol
    for i in range(len(a)):
        (sa, ea), (sb, eb) = a[i], b[i]
        if i == 0:
            if ea != eb or abs(sa - sb) > tol:
                return False
        elif i == len(a) - 1:
            if sa != sb or abs(ea - eb) > tol:
                return False
        elif (sa, ea) != (sb, eb):
            return False
    return True


def find_overlapping_gene(chrom, strand, start, end, by_loc):
    """Return gid of the gene on the same strand as the miniprot model that overlaps [start, end].
    Returns the highest-overlap gene if multiple. Strand must match."""
    best_gid = None
    best_overlap = 0
    intervals = by_loc.get((chrom, strand), [])
    starts = [iv[0] for iv in intervals]
    hi = bisect.bisect_right(starts, end)
    for i in range(hi):
        gs, ge, gid = intervals[i]
        if ge >= start:
            ov = min(end, ge) - max(start, gs) + 1
            if ov > best_overlap:
                best_overlap = ov
                best_gid = gid
    return best_gid


def main():
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} current.gff3 busco.gff missing_ids.txt output.gff3")

    current_path  = sys.argv[1]
    busco_path    = sys.argv[2]
    missing_path  = sys.argv[3]
    out_path      = sys.argv[4]

    missing_ids = load_missing(missing_path)
    sys.stderr.write(f"Missing BUSCOs to recover: {len(missing_ids)}\n")

    best_models = load_busco_models(busco_path, missing_ids)

    raw, gene_meta, by_loc, gene_order, gene_cds_sigs = load_annotation(current_path)
    gene_cds_sigs = defaultdict(set, {k: set(v) for k, v in gene_cds_sigs.items()})

    # For each missing BUSCO, decide where to attach it
    # gene_additions: gid → list of (mrna_gff_lines)  — isoforms to append after gene block
    # new_genes: list of (gene_line, [mrna_lines])  — standalone new genes
    gene_additions = defaultdict(list)   # gid → [lines to append]
    new_genes = []
    gene_expand = {}     # gid → (new_start, new_end)
    mrna_counter = [0]
    gene_counter = [0]

    n_isoform = n_new = n_notfound = n_dup_skipped = 0

    for bid, models in best_models.items():
        for model in models:   # add all models for each BUSCO (best-first)
            chrom, strand, start, end, score, orig_id, cds_rows = model

            gid = find_overlapping_gene(chrom, strand, start, end, by_loc)

            # Build mRNA + exon + CDS lines
            cds_rows_sorted = sorted(cds_rows, key=lambda c: int(c[3]))
            cds_list = [(int(c[3]), int(c[4])) for c in cds_rows_sorted]
            cds_sig = tuple(sorted(cds_list))

            if gid is not None and (
                cds_sig in gene_cds_sigs[gid]
                or any(_near_dup_cds(cds_list, sorted(sig)) for sig in gene_cds_sigs[gid])
            ):
                # An existing isoform of this gene already has the exact same
                # (or near-identical, modulo a partial start/stop codon) CDS —
                # the BUSCO is already represented in the protein set.
                # Attaching a duplicate would give table2asn two CDS features
                # at the same location, which antiSMASH rejects outright.
                n_dup_skipped += 1
                continue

            mrna_counter[0] += 1
            tid = f"busco_iso_{mrna_counter[0]:06d}"

            mrna_lines = []
            mrna_lines.append(
                f"{chrom}\tcompleasm\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t"
                f"ID={tid};alias=busco_{bid};Name=busco_{bid}_{orig_id}\n"
            )
            for i, c in enumerate(cds_rows_sorted):
                mrna_lines.append(
                    f"{chrom}\tcompleasm\texon\t{c[3]}\t{c[4]}\t.\t{strand}\t.\t"
                    f"ID={tid}.exon{i+1};Parent={tid}\n"
                )
            for i, c in enumerate(cds_rows_sorted):
                mrna_lines.append(
                    f"{chrom}\tcompleasm\tCDS\t{c[3]}\t{c[4]}\t.\t{strand}\t{c[7]}\t"
                    f"ID=cds.{tid}.{i+1};Parent={tid}\n"
                )

            if gid is not None:
                # Attach as isoform: update Parent, expand gene if needed
                mrna_lines[0] = mrna_lines[0].rstrip('\n').replace(
                    f"ID={tid};", f"ID={tid};Parent={gid};"
                ) + '\n'
                gene_additions[gid].append(mrna_lines)
                gene_cds_sigs[gid].add(cds_sig)
                # Track expansion needed
                cur = gene_meta[gid]
                new_s = min(cur['start'], start)
                new_e = max(cur['end'], end)
                if gid not in gene_expand:
                    gene_expand[gid] = [new_s, new_e]
                else:
                    gene_expand[gid][0] = min(gene_expand[gid][0], new_s)
                    gene_expand[gid][1] = max(gene_expand[gid][1], new_e)
                n_isoform += 1
            else:
                # No overlapping same-strand gene — add as standalone
                gene_counter[0] += 1
                new_gid = f"busco_gene_r{gene_counter[0]:06d}"
                gene_line = (
                    f"{chrom}\tcompleasm\tgene\t{start}\t{end}\t{score}\t{strand}\t.\t"
                    f"ID={new_gid};Name=busco_{bid}\n"
                )
                for j, l in enumerate(mrna_lines):
                    mrna_lines[j] = l.replace(f"ID={tid};", f"ID={tid};Parent={new_gid};")
                new_genes.append((gene_line, mrna_lines))
                # Add to by_loc so subsequent BUSCOs don't double-add
                by_loc.setdefault((chrom, strand), []).append((start, end, new_gid))
                by_loc[(chrom, strand)].sort()
                n_new += 1

    sys.stderr.write(f"Added as isoforms: {n_isoform}\n")
    sys.stderr.write(f"Added as new genes: {n_new}\n")
    sys.stderr.write(f"Skipped (CDS-duplicate of existing isoform): {n_dup_skipped}\n")

    # ── Write output ───────────────────────────────────────────────────────────
    # Track which genes we've closed (written all children)
    # Strategy: write raw, when we hit a gene line, check gene_expand and emit expansion.
    # After all children of that gene are written, append gene_additions.
    # After raw, append new_genes.

    cur_gid = None
    next_gene_idx = 0
    # Build gene-to-raw-end-line mapping: find line index of last child before next gene
    # Simpler: stream and track cur_gid; when gene changes, emit additions for previous gene

    with open(out_path, 'w') as out:
        cur_gid = None

        def close_gene():
            if cur_gid and cur_gid in gene_additions:
                for mrna_lines in gene_additions[cur_gid]:
                    for l in mrna_lines:
                        out.write(l)

        for line in raw:
            if line.startswith('#') or not line.strip():
                out.write(line)
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                out.write(line)
                continue
            feat = c[2]
            a = attr(c[8])

            if feat == 'gene':
                close_gene()
                gid = a.get('ID', '')
                cur_gid = gid
                if gid in gene_expand:
                    ns, ne = gene_expand[gid]
                    c[3] = str(ns)
                    c[4] = str(ne)
                    out.write('\t'.join(c) + '\n')
                else:
                    out.write(line)
            else:
                out.write(line)

        close_gene()   # last gene in file

        # Append new standalone genes
        if new_genes:
            out.write("##\n## --- Busco standalone new genes (recovery pass) ---\n##\n")
            for gene_line, mrna_lines in new_genes:
                out.write(gene_line)
                for l in mrna_lines:
                    out.write(l)

    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
