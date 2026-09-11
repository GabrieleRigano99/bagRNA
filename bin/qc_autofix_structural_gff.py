#!/usr/bin/env python3
"""QC gate for a structural annotation GFF3 — catches and auto-fixes the
same three NCBI table2asn error classes without needing a table2asn round
trip, so it can run right after any gene-set change (isoform addition,
BUSCO recovery, tRNA recovery, etc.) instead of requiring a manual
discovery cycle against the real validator.

Checks (all done locally against the genome FASTA):
  1. SEQ_INST.StopInProtein / SEQ_FEAT.InternalStop — translate every CDS,
     flag any with a stop codon before the final codon.
  2. SEQ_FEAT.ShortIntron — flag any CDS-internal intron < 10nt (the
     spliceosome's physical floor; anything smaller cannot be real).
  3. SEQ_FEAT.BadCDScomponentOverlapTRNA (and the rRNA/ncRNA equivalent) —
     flag any coding gene that overlaps a tRNA/rRNA/ncRNA AND has zero
     support across every functional-annotation column (requires
     --functional-tsv from merge_functional_annotations.py; skipped if
     not given).

Fix hierarchy (safest first, same logic validated in this project's
manual fix pass):
  a. If a flagged transcript has a clean sibling isoform, drop it —
     zero information loss, the gene is still represented.
  b. ShortIntron with no clean sibling: merge the two flanking CDS/exon
     segments (the "intron" bases become coding), recompute downstream
     phases, re-translate. Keep the merge if it's now clean.
  c. Whatever still has a real internal stop and no clean sibling: mark
     `pseudo=true` on gene/mRNA/CDS — the NCBI-sanctioned way to represent
     a genuine premature-stop CDS, rather than fabricating a "fixed"
     protein that might not reflect real biology.
  d. A coding gene fully overlapping a tRNA/rRNA/ncRNA with zero
     functional evidence anywhere: drop the coding gene, keep the ncRNA
     (tRNA/rRNA calls are far more reliable than an unsupported ORF call).

Scope note: short-intron detection/repair operates on CDS-internal
introns only (not UTR-internal), since that covers every case seen in
this project's annotation. UTR-internal micro-introns would need a
separate pass if they ever show up.

Usage:
    qc_autofix_structural_gff.py --gff struct.gff3 --fasta genome.fa \\
        [--functional-tsv functional_annotation.tsv] \\
        --out fixed.gff3 --report report.tsv
"""
import argparse
import sys
from collections import defaultdict

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
MIN_INTRON = 10
NC_TYPES = {'ncRNA', 'tRNA', 'rRNA', 'lncRNA'}
# Only tRNA/rRNA come from dedicated, high-confidence covariance-model/HMM
# tools (tRNAscan-SE, Barrnap) — they take precedence over an overlapping,
# unsupported coding-gene call. Plain 'ncRNA' here means Mikado-derived
# lncRNA calls, which are comparatively low-confidence and should NOT
# override a coding gene just because they overlap it (antisense lncRNA/
# coding-gene overlap is normal biology, not a sign the coding call is wrong).
HIGH_CONFIDENCE_NC_TYPES = {'tRNA', 'rRNA'}
CODING_TYPES = {'mRNA', 'transcript'}


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


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = ''.join(buf)
    return seqs


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    return ''.join(CODON_TABLE.get(seq[i:i+3].upper(), 'X') for i in range(0, len(seq) - 2, 3))


def cds_sequence(genome, chrom, segments, strand, start_phase=0):
    """segments: sorted (start, end) tuples, genomic order. start_phase is the
    GFF phase of whichever segment is translated FIRST (lowest coordinate on
    +, highest coordinate on -) — the number of bases to skip before the
    first complete codon. Needed for 5'-partial CDS (e.g. contig-edge genes),
    which otherwise translate out of frame and produce false internal stops.
    """
    seq = ''.join(genome[chrom][s - 1:e] for s, e in sorted(segments))
    if strand == '-':
        seq = revcomp(seq)
    return seq[start_phase:]


def first_segment_phase(cds_triples, strand):
    """cds_triples: (start, end, phase) tuples. Phase of whichever segment
    is translated first (lowest coord on +, highest coord on -)."""
    seg = max(cds_triples, key=lambda x: x[1]) if strand == '-' else min(cds_triples, key=lambda x: x[0])
    return seg[2]


def internal_stop_count(protein):
    body = protein[:-1] if protein.endswith('*') else protein
    return body.count('*')


# All-empty-evidence columns from merge_functional_annotations.py's TSV
EVIDENCE_COLS = [
    'product', 'GO_terms', 'EC_numbers', 'KEGG_KO', 'InterPro_accessions',
    'Pfam_domains', 'PANTHER', 'CAZyme_family', 'MEROPS_hit',
    'PHIbase_accession', 'BGC_cluster_type', 'SMCOG',
]


def load_functional_evidence(path):
    """gene_id -> True if it has zero support across every evidence column."""
    if not path:
        return {}
    no_evidence = {}
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) <= idx.get('gene_id', 0):
                continue
            gid = c[0]
            empty = True
            for col in EVIDENCE_COLS:
                i = idx.get(col)
                if i is None or i >= len(c):
                    continue
                v = c[i].strip()
                if v and v != 'hypothetical protein' and v != 'hypothetical lncRNA':
                    empty = False
                    break
            no_evidence[gid] = empty
    return no_evidence


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--functional-tsv', help='functional_annotation.tsv from merge_functional_annotations.py '
                                             '— enables the tRNA/rRNA overlap check')
    p.add_argument('--out', required=True)
    p.add_argument('--report', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    genome = load_fasta(args.fasta)
    no_evidence = load_functional_evidence(args.functional_tsv)

    genes = []
    cur_gene = None
    with open(args.gff) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            feat = c[2]
            a = attr(c[8])
            if feat == 'gene':
                cur_gene = {'cols': c[:], 'attrs': a, 'transcripts': []}
                genes.append(cur_gene)
            elif feat in CODING_TYPES | NC_TYPES:
                if cur_gene is None:
                    continue
                t = {'old_id': a.get('ID', ''), 'type': feat, 'cols': c[:], 'attrs': a, 'children': []}
                cur_gene['transcripts'].append(t)
            else:
                if cur_gene is None or not cur_gene['transcripts']:
                    continue
                parent = a.get('Parent', '')
                for t in reversed(cur_gene['transcripts']):
                    if t['old_id'] == parent:
                        t['children'].append({'type': feat, 'cols': c[:], 'attrs': a})
                        break

    gene_mrnas = defaultdict(list)
    gene_is_coding = {}
    for g in genes:
        gid = g['attrs'].get('ID', '')
        types = {t['type'] for t in g['transcripts']}
        gene_is_coding[gid] = bool(types & CODING_TYPES)
        for t in g['transcripts']:
            if t['type'] in CODING_TYPES:
                gene_mrnas[gid].append(t['old_id'])

    report = []

    # ── 1+2: local detection — internal stops and short introns ────────────
    bad_stop = set()
    short_intron = {}   # tid -> (gap_start, gap_end)
    for g in genes:
        for t in g['transcripts']:
            if t['type'] not in CODING_TYPES:
                continue
            cds3 = sorted((int(ch['cols'][3]), int(ch['cols'][4]), int(ch['cols'][7]))
                          for ch in t['children'] if ch['type'] == 'CDS')
            if not cds3:
                continue
            cds = [(s, e) for s, e, _ in cds3]
            chrom, strand = t['cols'][0], t['cols'][6]
            phase0 = first_segment_phase(cds3, strand)
            protein = translate(cds_sequence(genome, chrom, cds, strand, phase0))
            if internal_stop_count(protein) > 0:
                bad_stop.add(t['old_id'])
            for i in range(len(cds) - 1):
                gap = cds[i + 1][0] - cds[i][1] - 1
                if 0 < gap < MIN_INTRON:
                    short_intron[t['old_id']] = (cds[i][1], cds[i + 1][0])
                    break   # one flagged gap per transcript is enough to drive the same fix as before

    # ── ShortIntron fix attempt ──────────────────────────────────────────────
    fixed_clean = set()
    for g in genes:
        for t in g['transcripts']:
            tid = t['old_id']
            if tid not in short_intron:
                continue
            gap_s, gap_e = short_intron[tid]
            cds = sorted([(int(ch['cols'][3]), int(ch['cols'][4]), ch) for ch in t['children'] if ch['type'] == 'CDS'])
            phase_lut = {(int(ch['cols'][3]), int(ch['cols'][4])): int(ch['cols'][7])
                         for ch in t['children'] if ch['type'] == 'CDS'}
            exons = sorted([(int(ch['cols'][3]), int(ch['cols'][4]), ch) for ch in t['children'] if ch['type'] == 'exon'])
            m = next((i for i in range(len(cds) - 1) if cds[i][1] == gap_s and cds[i + 1][0] == gap_e), None)
            if m is None:
                report.append((tid, 'ShortIntron', 'ERROR: could not locate flanking CDS segments'))
                continue
            i, j = m, m + 1
            new_cds = [(s, e) for k, (s, e, _) in enumerate(cds) if k not in (i, j)] + [(cds[i][0], cds[j][1])]
            new_cds.sort()
            em = next((k for k in range(len(exons) - 1) if exons[k][1] == gap_s and exons[k + 1][0] == gap_e), None)
            new_exon = [(s, e) for k, (s, e, _) in enumerate(exons)]
            if em is not None:
                new_exon = [(s, e) for k, (s, e) in enumerate(new_exon) if k not in (em, em + 1)] + [(exons[em][0], exons[em + 1][1])]
                new_exon.sort()

            chrom, strand = t['cols'][0], t['cols'][6]
            # Anchor coordinate of the translation-first NEW segment maps back to
            # whichever original segment supplied that end — its phase carries over
            # (merging only extends the OTHER end of a segment, never the anchor end).
            anchor = max(new_cds, key=lambda x: x[1])[1] if strand == '-' else min(new_cds, key=lambda x: x[0])[0]
            phase0 = next((p for (s, e), p in phase_lut.items() if (strand == '-' and e == anchor) or (strand == '+' and s == anchor)), 0)
            protein = translate(cds_sequence(genome, chrom, new_cds, strand, phase0))
            n_stops = internal_stop_count(protein)
            t['_new_cds'], t['_new_exon'] = new_cds, new_exon
            if n_stops == 0:
                fixed_clean.add(tid)
                report.append((tid, 'ShortIntron', f'merged flanking exons ({gap_e - gap_s + 1}bp gap) -> clean ORF'))
            else:
                report.append((tid, 'ShortIntron', f'merged flanking exons but {n_stops} internal stop(s) remain'))

    still_bad = (bad_stop | set(short_intron)) - fixed_clean

    def has_clean_sibling(gid, tid):
        return any(sib != tid and sib not in still_bad for sib in gene_mrnas.get(gid, []))

    drop_transcripts, pseudogene_transcripts = set(), set()
    for g in genes:
        gid = g['attrs'].get('ID', '')
        for t in g['transcripts']:
            tid = t['old_id']
            if tid not in still_bad:
                continue
            cause = 'StopInProtein' if tid in bad_stop else 'ShortIntron (merge still broken)'
            if has_clean_sibling(gid, tid):
                drop_transcripts.add(tid)
                report.append((tid, cause, 'dropped — clean sibling isoform exists'))
            else:
                pseudogene_transcripts.add(tid)
                report.append((tid, cause, 'marked pseudogene — sole isoform, no clean alternative'))

    # ── 3: tRNA/rRNA vs zero-evidence coding gene overlap ────────────────────
    # Deliberately excludes plain ncRNA/lncRNA — see HIGH_CONFIDENCE_NC_TYPES.
    drop_genes = set()
    if no_evidence:
        nc_spans = []   # (chrom, start, end, gid)
        for g in genes:
            gid = g['attrs'].get('ID', '')
            if any(t['type'] in HIGH_CONFIDENCE_NC_TYPES for t in g['transcripts']):
                nc_spans.append((g['cols'][0], int(g['cols'][3]), int(g['cols'][4]), gid))
        for g in genes:
            gid = g['attrs'].get('ID', '')
            if not gene_is_coding.get(gid):
                continue
            if not no_evidence.get(gid, False):
                continue
            chrom, gs, ge = g['cols'][0], int(g['cols'][3]), int(g['cols'][4])
            for nchrom, ns, ne, ngid in nc_spans:
                if nchrom == chrom and gs <= ne and ns <= ge:
                    drop_genes.add(gid)
                    report.append((gid, 'NC-overlap', f'dropped — zero functional evidence, overlaps {ngid}'))
                    break

    # ── Write output ─────────────────────────────────────────────────────────
    def write_feature(out, cols, a):
        c = cols[:]
        c[8] = attr_str(a) + ';'
        out.write('\t'.join(c) + '\n')

    n_dropped_t = n_pseudo = n_dropped_g = 0
    with open(args.out, 'w') as out:
        out.write('##gff-version 3\n')
        for g in genes:
            gid = g['attrs'].get('ID', '')
            if gid in drop_genes:
                n_dropped_g += 1
                continue
            kept = [t for t in g['transcripts'] if t['old_id'] not in drop_transcripts]
            if not kept:
                continue
            gene_is_pseudo = any(t['old_id'] in pseudogene_transcripts for t in kept)
            g_attrs = dict(g['attrs'])
            if gene_is_pseudo:
                g_attrs['pseudo'] = 'true'
                g_attrs['pseudogene'] = 'unknown'
                n_pseudo += 1
            write_feature(out, g['cols'], g_attrs)

            for t in kept:
                tid = t['old_id']
                is_pseudo_here = tid in pseudogene_transcripts
                t_attrs = dict(t['attrs'])
                if is_pseudo_here:
                    t_attrs['pseudo'] = 'true'
                write_feature(out, t['cols'], t_attrs)

                if '_new_cds' in t:
                    for s, e in sorted(t['_new_exon']):
                        proto = next(ch for ch in t['children'] if ch['type'] == 'exon')
                        ec = proto['cols'][:]
                        ec[3], ec[4] = str(s), str(e)
                        write_feature(out, ec, dict(proto['attrs']))
                    ordered = sorted(t['_new_cds']) if t['cols'][6] == '+' else sorted(t['_new_cds'], reverse=True)
                    cum = 0
                    proto = next(ch for ch in t['children'] if ch['type'] == 'CDS')
                    for s, e in ordered:
                        phase = (3 - (cum % 3)) % 3
                        cc = proto['cols'][:]
                        cc[3], cc[4], cc[7] = str(s), str(e), str(phase)
                        ca = dict(proto['attrs'])
                        if is_pseudo_here:
                            ca['pseudo'] = 'true'
                        write_feature(out, cc, ca)
                        cum += e - s + 1
                else:
                    for ch in t['children']:
                        ca = dict(ch['attrs'])
                        if is_pseudo_here:
                            ca['pseudo'] = 'true'
                        write_feature(out, ch['cols'], ca)

            for t in g['transcripts']:
                if t['old_id'] in drop_transcripts:
                    n_dropped_t += 1

    with open(args.report, 'w') as rf:
        rf.write('id\tissue\taction\n')
        for tid, issue, action in report:
            rf.write(f"{tid}\t{issue}\t{action}\n")

    sys.stderr.write(f"Internal-stop transcripts found: {len(bad_stop)}\n")
    sys.stderr.write(f"Short-intron transcripts found: {len(short_intron)}\n")
    sys.stderr.write(f"  -> merged to a clean ORF: {len(fixed_clean)}\n")
    sys.stderr.write(f"Dropped transcripts (clean sibling existed): {n_dropped_t}\n")
    sys.stderr.write(f"Marked pseudogene (sole isoform, unfixable): {n_pseudo}\n")
    sys.stderr.write(f"Dropped zero-evidence genes overlapping ncRNA/tRNA/rRNA: {n_dropped_g}\n")
    sys.stderr.write(f"Wrote {args.out}\n")
    sys.stderr.write(f"Wrote report: {args.report}\n")


if __name__ == '__main__':
    main()
