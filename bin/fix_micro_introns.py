#!/usr/bin/env python3
"""Merge spurious sub-10bp "introns" into the adjacent CDS/exon.

Real spliceosomal introns are never this short (canonical minimum is
~30-70bp); a gap below NCBI's hard minimum (10 nt) between two CDS
segments of the same transcript is an ab initio boundary-calling
artifact (observed from ANNEVO primaries Mikado picked as-is), not a
genuine splice site. table2asn rejects the whole gene for it, so
recovering it here — before GFF_CLEAN_FILTER's NCBI QC — turns a
correctly-formed but falsely-split gene back into a kept one.

Safety: retaining a gap as coding only preserves the reading frame for
everything downstream if the gap length is a multiple of 3 (checked).
The only *new* codons created by the merge are the ones straddling the
old segment boundaries — any codon fully inside the original CDS is
unaffected and re-translates identically to before. So we only need to
validate that straddling region (leftover tail of the upstream segment
+ the whole gap + however many bases of the downstream segment are
needed to complete the last new codon) for an introduced stop codon.
If a gap fails either check, it is left untouched as a real intron —
this script can only rescue genes, never break a working one.

A merge can incidentally make one isoform's CDS become byte-identical to
a sibling isoform's (e.g. a protected Mikado primary with a spurious
micro-intron ends up matching a TransDecoder-derived isoform of the same
gene that never had the gap) — upstream dedup (clean_isoforms_gff.py)
ran before this merge existed, so it never saw the collision. table2asn
happily emits two CDS features at the exact same location, which
antiSMASH then rejects outright ("Multiple CDS features have the same
location"). So after merging, every gene's surviving isoforms are
re-checked for CDS-identical duplicates (independent of whether a merge
caused them) and only the one with the most UTR support (longest mRNA
span) is kept.

Usage:
    fix_micro_introns.py input.gff3 genome.fasta output.gff3 [--max-gap N]
"""
import sys
from collections import defaultdict

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def load_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def attr(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


def has_stop(seq):
    return any(seq[i:i + 3] in STOP_CODONS for i in range(0, len(seq) - 2, 3))


def seg_seq(genome, chrom, s, e, strand):
    # .upper(): genome FASTAs are soft-masked (lowercase repeats) — stop-codon
    # detection must not miss a masked "taa"/"tag"/"tga"
    seq = genome[chrom][s - 1:e].upper()
    return revcomp(seq) if strand == '-' else seq


def main():
    if len(sys.argv) < 4:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 genome.fasta output.gff3 [--max-gap N]")

    in_path, fasta_path, out_path = sys.argv[1:4]
    max_gap = 10
    if '--max-gap' in sys.argv:
        max_gap = int(sys.argv[sys.argv.index('--max-gap') + 1])

    genome = load_fasta(fasta_path)

    with open(in_path) as f:
        raw = f.readlines()

    # ── Collect CDS segments per mRNA ─────────────────────────────────────────
    cds_by_tx = defaultdict(list)
    meta = {}
    gene_of = {}
    span_of = {}
    cur_tid = None

    for line in raw:
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        feat = c[2]
        a = attr(c[8])
        if feat in ('mRNA', 'transcript'):
            cur_tid = a.get('ID', '')
            meta[cur_tid] = (c[0], c[6])
            gene_of[cur_tid] = a.get('Parent', '')
            span_of[cur_tid] = (int(c[3]), int(c[4]))
        elif feat == 'CDS':
            parent = a.get('Parent', cur_tid)
            cds_by_tx[parent].append((int(c[3]), int(c[4])))
        elif feat == 'gene':
            cur_tid = None

    # ── Decide which gaps are safe to merge ───────────────────────────────────
    merges = defaultdict(set)   # tid -> {i, ...} meaning order[i] merges with order[i+1]
    n_merged_genes = set()
    n_merged_gaps = 0

    for tid, segs in cds_by_tx.items():
        if len(segs) < 2:
            continue
        chrom, strand = meta.get(tid, (None, None))
        if chrom not in genome:
            continue
        segs_genomic = sorted(segs)
        order = segs_genomic if strand == '+' else list(reversed(segs_genomic))

        cum_len = 0  # CDS bases consumed strictly before the current segment
        for i in range(len(order) - 1):
            cur_s, cur_e = order[i]
            nxt_s, nxt_e = order[i + 1]
            cur_len = cur_e - cur_s + 1
            phase_in = cum_len % 3           # phase entering this segment
            cum_len += cur_len

            if strand == '+':
                gap_start, gap_end = cur_e + 1, nxt_s - 1
            else:
                gap_start, gap_end = nxt_e + 1, cur_s - 1
            gap_len = gap_end - gap_start + 1

            if not (0 <= gap_len < max_gap) or gap_len % 3 != 0:
                continue  # not a candidate, or unsafe (would shift downstream frame)

            r1 = (cur_len - phase_in) % 3    # leftover un-codon-completed tail of cur
            r2 = (3 - r1) % 3                # bases needed from nxt to complete that codon

            gap_seq = seg_seq(genome, chrom, gap_start, gap_end, strand)
            tail = seg_seq(genome, chrom, cur_e - r1 + 1, cur_e, strand) if r1 else ''
            head = seg_seq(genome, chrom, nxt_s, nxt_s + r2 - 1, strand) if r2 else ''
            straddle = tail + gap_seq + head

            if len(straddle) % 3 != 0 or has_stop(straddle):
                continue

            merges[tid].add(i)
            n_merged_genes.add(tid)
            n_merged_gaps += 1

    sys.stderr.write(
        f"Merging {n_merged_gaps} sub-{max_gap}bp gap(s) across {len(n_merged_genes)} "
        f"transcript(s) — frame-safe, stop-free retained-sequence merges only\n"
    )

    # ── Compute merged CDS coordinates per transcript (genomic order out) ─────
    new_segs = {}
    for tid in merges:
        chrom, strand = meta[tid]
        segs_genomic = sorted(cds_by_tx[tid])
        order = segs_genomic if strand == '+' else list(reversed(segs_genomic))

        merged = []
        i = 0
        while i < len(order):
            s, e = order[i]
            while i in merges[tid]:
                i += 1
                ns, ne = order[i]
                s, e = min(s, ns), max(e, ne)
            merged.append((s, e))
            i += 1
        new_segs[tid] = sorted(merged)

    # ── Drop isoforms whose (possibly just-merged) CDS now duplicates a
    #    sibling's, keeping the one with the most UTR support ────────────────
    final_cds = {
        tid: tuple(new_segs[tid]) if tid in new_segs else tuple(sorted(segs))
        for tid, segs in cds_by_tx.items()
    }

    groups = defaultdict(list)
    for tid, cds_key in final_cds.items():
        groups[(gene_of.get(tid), cds_key)].append(tid)

    drop_tid = set()
    for (gid, cds_key), tids in groups.items():
        if gid is None or len(tids) < 2:
            continue
        tids.sort(key=lambda t: span_of[t][0] - span_of[t][1])  # widest span first
        drop_tid.update(tids[1:])

    if drop_tid:
        sys.stderr.write(
            f"Dropping {len(drop_tid)} isoform(s) whose CDS now duplicates a "
            f"sibling's after micro-intron merging (kept the longest-span one)\n"
        )

    # ── Rewrite the GFF: replace CDS/exon lines for merged transcripts,
    #    drop transcripts (and all their children) that became CDS-duplicates ──
    with open(out_path, 'w') as out:
        cur_tid = None
        skip_children_of = None
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
                cur_tid = None
                skip_children_of = None
                out.write(line)
            elif feat in ('mRNA', 'transcript'):
                cur_tid = a.get('ID', '')
                if cur_tid in drop_tid:
                    skip_children_of = cur_tid
                    continue
                out.write(line)
                if cur_tid in new_segs:
                    skip_children_of = cur_tid
                    chrom, strand = meta[cur_tid]
                    segs = new_segs[cur_tid]
                    for j, (s, e) in enumerate(segs, 1):
                        out.write(f"{chrom}\t{c[1]}\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                                  f"ID={cur_tid}.mifix_exon{j};Parent={cur_tid}\n")
                    seg_order = segs if strand == '+' else list(reversed(segs))
                    phase = 0
                    for j, (s, e) in enumerate(seg_order, 1):
                        out.write(f"{chrom}\t{c[1]}\tCDS\t{s}\t{e}\t.\t{strand}\t{phase}\t"
                                  f"ID={cur_tid}.mifix_cds{j};Parent={cur_tid}\n")
                        phase = (3 - ((e - s + 1 - phase) % 3)) % 3
                else:
                    skip_children_of = None
            elif skip_children_of and a.get('Parent') == skip_children_of:
                continue  # original pre-merge CDS/exon, or any child of a dropped duplicate
            else:
                out.write(line)

    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
