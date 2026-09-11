#!/usr/bin/env python3
"""Clean a Mikado-derived GFF3 after isoform addition.

Four operations (applied in order, iterated until convergence):
  1. Remove isoforms with a duplicate exon structure within the same gene.
  2. Remove isoforms that are redundant within the same gene:
     a. Same splice junctions as another isoform → drop the shorter.
     b. All exons contained within another isoform → drop the contained one.
  3. Update gene feature coordinates to encompass all kept child mRNAs.
  4. Resolve all same-strand gene overlaps:
     - UTR-only overlap: clip the UTR(s) so the genes no longer touch.
     - CDS-level overlap: drop the lower-scoring gene (Mikado score col).

Usage:
    clean_isoforms_gff.py input.gff3 output.gff3
"""
import sys
from collections import defaultdict


def attr(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


def _junctions(exons):
    s = sorted(exons)
    return frozenset((s[i][1], s[i + 1][0]) for i in range(len(s) - 1))


def _contained_in(a, b):
    """True if every exon of 'a' is spatially contained within some exon of 'b'
    and the span of 'a' is within the span of 'b'."""
    if a['start'] < b['start'] or a['end'] > b['end']:
        return False
    b_exons = sorted(b['exons'])
    for sa, ea in a['exons']:
        if not any(sb <= sa and ea <= eb for sb, eb in b_exons):
            return False
    return True


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 output.gff3")

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path) as f:
        raw = f.readlines()

    # ── First pass: collect structure ─────────────────────────────────────────
    gene_lines = {}    # gid → line index
    gene_meta  = {}    # gid → {chrom, strand, score}
    mrna_lines = {}    # tid → line index
    mrnas      = {}    # tid → {gene, start, end, exons, cds}
    cur_tid    = None

    for i, line in enumerate(raw):
        if line.startswith('#') or not line.strip():
            continue
        c = line.split('\t')
        if len(c) < 9:
            continue
        a    = attr(c[8])
        feat = c[2]

        if feat == 'gene':
            gid = a.get('ID', '')
            gene_lines[gid] = i
            try:
                score = float(c[5])
            except ValueError:
                score = 0.0
            gene_meta[gid] = {'chrom': c[0], 'strand': c[6], 'score': score}
            cur_tid = None

        elif feat in ('mRNA', 'ncRNA', 'lncRNA', 'transcript'):
            parent = a.get('Parent', '')
            if parent not in gene_meta:
                cur_tid = None   # orphan transcript — leave as-is
                continue
            cur_tid = a.get('ID', '')
            mrna_lines[cur_tid] = i
            mrnas[cur_tid] = {
                'gene':  parent,
                'start': int(c[3]), 'end': int(c[4]),
                'exons': [], 'cds': [],
                'source': c[1],   # track source for protection
            }

        elif feat == 'exon':
            parent = a.get('Parent', cur_tid)
            if parent in mrnas:
                mrnas[parent]['exons'].append((int(c[3]), int(c[4])))

        elif feat == 'CDS':
            parent = a.get('Parent', cur_tid)
            if parent in mrnas:
                mrnas[parent]['cds'].append((int(c[3]), int(c[4])))

    # ── Deduplicate isoforms by exon structure ────────────────────────────────
    # "Protected" = original annotation mRNAs (source != transdecoder / compleasm).
    # These are never dropped regardless of duplication or containment.
    ADDED_SOURCES = {'transdecoder', 'compleasm'}
    protected_tids = {tid for tid, m in mrnas.items()
                      if m.get('source', '') not in ADDED_SOURCES}

    by_gene   = defaultdict(list)
    for tid, m in mrnas.items():
        by_gene[m['gene']].append(tid)

    drop_mrna = set()
    for gid, tids in by_gene.items():
        seen = {}    # exon_key → tid that claimed it
        for tid in tids:
            key = tuple(sorted(mrnas[tid]['exons']))
            if key in seen:
                # Keep protected over non-protected; otherwise keep first seen
                prev = seen[key]
                if tid in protected_tids and prev not in protected_tids:
                    drop_mrna.add(prev)
                    seen[key] = tid   # replace with protected
                else:
                    drop_mrna.add(tid)
            else:
                seen[key] = tid

    sys.stderr.write(f"Dropping {len(drop_mrna)} duplicate-exon-structure isoforms\n")

    # ── Drop redundant isoforms: same junctions or contained in another ───────
    n_redundant = 0
    for gid, tids in by_gene.items():
        kept = [t for t in tids if t not in drop_mrna]
        if len(kept) < 2:
            continue
        # Sort: protected first, then by longest span.
        # This ensures protected mRNAs are always 'ta' (the keeper) when there is a tie.
        kept.sort(key=lambda t: (t not in protected_tids,
                                 -(mrnas[t]['end'] - mrnas[t]['start'])))
        for i, ta in enumerate(kept):
            if ta in drop_mrna:
                continue
            ja = _junctions(mrnas[ta]['exons'])
            for tb in kept[i + 1:]:
                if tb in drop_mrna:
                    continue
                jb = _junctions(mrnas[tb]['exons'])
                # Same splice junctions (multi-exonic only) → drop shorter (tb).
                # Never drop a protected mRNA.
                if ja and ja == jb and tb not in protected_tids:
                    drop_mrna.add(tb)
                    n_redundant += 1
                # tb's exons entirely contained in ta → drop tb (if not protected)
                elif _contained_in(mrnas[tb], mrnas[ta]) and tb not in protected_tids:
                    drop_mrna.add(tb)
                    n_redundant += 1

    sys.stderr.write(f"Dropping {n_redundant} redundant isoforms (same junctions or contained)\n")

    # ── Drop isoforms with an identical CDS to another isoform of the same gene ──
    # Isoforms can differ only in UTR extent (different exon/junction sets from
    # the passes above) while encoding the exact same protein. table2asn emits
    # one CDS feature per surviving mRNA, so duplicate CDS spans within a gene
    # produce multiple CDS features at the same genomic location — antiSMASH
    # rejects this ("Multiple CDS features have the same location"). Since the
    # duplicates are protein-identical, only the isoform with the most UTR
    # coverage (longest mRNA span) needs to survive; ties prefer protected.
    n_cds_dup = 0
    for gid, tids in by_gene.items():
        kept = [t for t in tids if t not in drop_mrna]
        if len(kept) < 2:
            continue
        cds_groups = defaultdict(list)
        for t in kept:
            if not mrnas[t]['cds']:
                continue   # noncoding transcript — not part of CDS dedup
            cds_groups[tuple(sorted(mrnas[t]['cds']))].append(t)
        for cds_key, group in cds_groups.items():
            if len(group) < 2:
                continue
            group.sort(key=lambda t: (t not in protected_tids,
                                      -(mrnas[t]['end'] - mrnas[t]['start'])))
            for t in group[1:]:
                drop_mrna.add(t)
                n_cds_dup += 1

    sys.stderr.write(f"Dropping {n_cds_dup} CDS-duplicate isoforms (identical protein, redundant UTR variant)\n")

    # ── Expand gene coordinates to encompass all kept mRNAs ──────────────────
    gene_new_coords = {}   # gid → [start, end]  (mutable)
    for gid, tids in by_gene.items():
        kept = [t for t in tids if t not in drop_mrna]
        if not kept:
            continue
        gene_new_coords[gid] = [
            min(mrnas[t]['start'] for t in kept),
            max(mrnas[t]['end']   for t in kept),
        ]

    expanded = sum(
        1 for gid, (ns, ne) in gene_new_coords.items()
        if (int(raw[gene_lines[gid]].split('\t')[3]) != ns or
            int(raw[gene_lines[gid]].split('\t')[4]) != ne)
    )
    sys.stderr.write(f"Expanding coordinates for {expanded} genes\n")

    # ── Resolve same-strand gene overlaps (iterates until convergence) ────────
    gene_is_coding = {
        gid: any(bool(mrnas[t]['cds']) for t in tids)
        for gid, tids in by_gene.items()
    }

    locus_genes = defaultdict(list)   # (chrom, strand) → [gid, ...]
    for gid in gene_new_coords:
        m = gene_meta[gid]
        locus_genes[(m['chrom'], m['strand'])].append(gid)

    drop_gene = set()
    mrna_clip = {}   # tid → [clip_start, clip_end]
    n_clipped  = 0
    n_dropped  = 0
    n_nc_drops = 0

    changed = True
    while changed:
        changed = False
        for gids in locus_genes.values():
            # Build sorted list of currently valid genes
            valid = sorted(
                [gene_new_coords[g] + [g] for g in gids if g not in drop_gene]
            )   # each entry: [start, end, gid]

            i = 0
            while i < len(valid) - 1:
                s1, e1, gid1 = valid[i]
                s2, e2, gid2 = valid[i + 1]

                if s2 > e1:          # no overlap
                    i += 1
                    continue

                # Coding vs noncoding: always drop the noncoding gene
                coding1 = gene_is_coding.get(gid1, False)
                coding2 = gene_is_coding.get(gid2, False)
                if coding1 and not coding2:
                    drop_gene.add(gid2)
                    valid.pop(i + 1)
                    n_nc_drops += 1
                    changed = True
                    continue
                if coding2 and not coding1:
                    drop_gene.add(gid1)
                    valid.pop(i)
                    n_nc_drops += 1
                    changed = True
                    continue

                kept1 = [t for t in by_gene.get(gid1, []) if t not in drop_mrna]
                kept2 = [t for t in by_gene.get(gid2, []) if t not in drop_mrna]

                cds_end1 = max(
                    (max((e for _, e in mrnas[t]['cds']), default=mrnas[t]['end'])
                     for t in kept1),
                    default=e1
                )
                cds_start2 = min(
                    (min((s for s, _ in mrnas[t]['cds']), default=mrnas[t]['start'])
                     for t in kept2),
                    default=s2
                )

                cut_right = s2 - 1   # clip gid1 here
                cut_left  = e1 + 1   # clip gid2 here

                if cds_end1 <= cut_right:
                    for t in kept1:
                        if mrnas[t]['end'] > cut_right:
                            clip = mrna_clip.setdefault(t, [mrnas[t]['start'], mrnas[t]['end']])
                            clip[1] = min(clip[1], cut_right)
                    new_e1 = max(min(mrnas[t]['end'], cut_right) for t in kept1)
                    gene_new_coords[gid1][1] = new_e1
                    valid[i][1] = new_e1
                    n_clipped += 1
                    changed = True
                    i += 1

                elif cds_start2 > e1:
                    for t in kept2:
                        if mrnas[t]['start'] < cut_left:
                            clip = mrna_clip.setdefault(t, [mrnas[t]['start'], mrnas[t]['end']])
                            clip[0] = max(clip[0], cut_left)
                    new_s2 = min(max(mrnas[t]['start'], cut_left) for t in kept2)
                    gene_new_coords[gid2][0] = new_s2
                    valid[i + 1][0] = new_s2
                    n_clipped += 1
                    changed = True
                    i += 1

                else:
                    # CDS-level overlap — drop lower-scoring gene
                    score1 = gene_meta[gid1]['score']
                    score2 = gene_meta[gid2]['score']
                    if score1 >= score2:
                        drop_gene.add(gid2)
                        valid.pop(i + 1)
                    else:
                        drop_gene.add(gid1)
                        valid.pop(i)
                    n_dropped += 1
                    changed = True
                    # don't increment i — re-check same position

    sys.stderr.write(f"Resolved {n_clipped} overlapping pairs by UTR clipping\n")
    sys.stderr.write(f"Dropped {n_nc_drops} noncoding genes overlapping coding genes\n")
    sys.stderr.write(f"Dropped {n_dropped} genes with CDS-level overlap (kept higher score)\n")

    # ── Second pass: write corrected GFF ─────────────────────────────────────
    skip_gene   = False   # True: we are inside a dropped gene — skip everything
    skip_parent = None    # mRNA ID whose children to skip

    with open(out_path, 'w') as out:
        for line in raw:
            if line.startswith('#') or not line.strip():
                if not skip_gene:
                    out.write(line)
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                if not skip_gene:
                    out.write(line)
                continue
            feat = c[2]
            a    = attr(c[8])

            if feat == 'gene':
                gid = a.get('ID', '')
                if gid in drop_gene:
                    skip_gene = True
                    skip_parent = None
                    continue
                skip_gene   = False
                skip_parent = None
                if gid in gene_new_coords:
                    ns, ne = gene_new_coords[gid]
                    c[3], c[4] = str(ns), str(ne)
                out.write('\t'.join(c) + '\n')

            elif skip_gene:
                continue

            elif feat in ('mRNA', 'ncRNA', 'lncRNA', 'transcript'):
                tid = a.get('ID', '')
                if tid in drop_mrna:
                    skip_parent = tid
                    continue
                skip_parent = None
                if tid in mrna_clip:
                    cs, ce = mrna_clip[tid]
                    c[3] = str(max(int(c[3]), cs))
                    c[4] = str(min(int(c[4]), ce))
                out.write('\t'.join(c) + '\n')

            else:
                parent = a.get('Parent', '')
                if parent == skip_parent:
                    continue
                if parent in mrna_clip:
                    cs, ce = mrna_clip[parent]
                    fs, fe = int(c[3]), int(c[4])
                    ns, ne = max(fs, cs), min(fe, ce)
                    if ns > ne:
                        continue   # entirely clipped away
                    if ns != fs or ne != fe:
                        c[3], c[4] = str(ns), str(ne)
                        out.write('\t'.join(c) + '\n')
                        continue
                out.write(line)


if __name__ == '__main__':
    main()
