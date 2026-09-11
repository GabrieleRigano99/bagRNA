#!/usr/bin/env python3
"""Fix the gene models table2asn flagged as StopInProtein/InternalStop/ShortIntron.

Strategy (safest-first):
  1. If the flagged transcript has a clean sibling isoform on the same gene,
     drop the flagged one. Zero information loss — the gene is still
     represented correctly.
  2. ShortIntron, no clean sibling: these gaps are 6-9 nt — physically too
     small for the spliceosome to excise (NCBI's own 10nt floor). Merge the
     two flanking CDS/exon segments into one continuous segment (the
     "intron" bases become coding), recompute downstream CDS phases, and
     re-translate to confirm the fix didn't introduce a new stop.
  3. Whatever is still broken after (1) and (2) — a transcript with a real
     internal stop and no clean sibling — gets marked `pseudo=true` on its
     gene/mRNA/CDS. This is the standard NCBI-sanctioned way to represent a
     CDS with a genuine premature stop, instead of fabricating a "repaired"
     protein sequence that might not reflect real biology.

Usage:
    fix_flagged_genes.py --gff struct_final_v3_renamed.gff3 --fasta genome.fa \\
        --stopinprotein bad_stopinprotein.txt --shortintron shortintron_coords.tsv \\
        --out fixed.gff3 --report report.tsv
"""
import argparse
import re
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


def translate(cds_seq):
    prot = []
    for i in range(0, len(cds_seq) - 2, 3):
        codon = cds_seq[i:i+3].upper()
        prot.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(prot)


def cds_sequence(genome, chrom, segments, strand):
    """segments: sorted list of (start, end), 1-based inclusive, genomic order."""
    parts = [genome[chrom][s - 1:e] for s, e in segments]
    seq = ''.join(parts)
    return revcomp(seq) if strand == '-' else seq


def internal_stop_count(protein):
    body = protein[:-1] if protein.endswith('*') else protein
    return body.count('*')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--stopinprotein', required=True, help='one transcript ID per line')
    p.add_argument('--shortintron', required=True, help='TSV: transcript_id\\tgap_start\\tgap_end')
    p.add_argument('--out', required=True)
    p.add_argument('--report', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    genome = load_fasta(args.fasta)
    bad_stop = set(open(args.stopinprotein).read().split())
    short_intron = {}   # tid -> (gap_start, gap_end)
    for line in open(args.shortintron):
        c = line.split()
        if len(c) >= 3:
            short_intron[c[0]] = (int(c[1]), int(c[2]))

    # ── Parse GFF into gene -> transcripts -> children ──────────────────────
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
            elif feat in ('mRNA', 'ncRNA', 'tRNA', 'rRNA', 'transcript'):
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
    for g in genes:
        gid = g['attrs'].get('ID', '')
        for t in g['transcripts']:
            if t['type'] in ('mRNA', 'transcript'):
                gene_mrnas[gid].append(t['old_id'])

    report = []

    # ── Pass 1: try merging short introns, check if translation is clean ────
    fixed_clean = set()     # transcripts where the merge produced a clean ORF
    fixed_still_bad = set()  # transcripts where merge happened but a stop remains

    for g in genes:
        for t in g['transcripts']:
            tid = t['old_id']
            if tid not in short_intron:
                continue
            gap_s, gap_e = short_intron[tid]
            cds = sorted(
                [(int(ch['cols'][3]), int(ch['cols'][4]), ch) for ch in t['children'] if ch['type'] == 'CDS'],
                key=lambda x: x[0]
            )
            exons = sorted(
                [(int(ch['cols'][3]), int(ch['cols'][4]), ch) for ch in t['children'] if ch['type'] == 'exon'],
                key=lambda x: x[0]
            )
            merged_cds = None
            for i in range(len(cds) - 1):
                if cds[i][1] == gap_s and cds[i + 1][0] == gap_e:
                    merged_cds = (i, i + 1)
                    break
            if merged_cds is None:
                report.append((tid, 'ShortIntron', 'ERROR: could not locate flanking CDS segments'))
                continue
            i, j = merged_cds
            new_start, new_end = cds[i][0], cds[j][1]
            new_cds_segments = [(s, e) for k, (s, e, _) in enumerate(cds) if k not in (i, j)]
            new_cds_segments.append((new_start, new_end))
            new_cds_segments.sort()

            # mirror the merge on the exon list (same coordinates at this junction)
            merged_exon = None
            for k in range(len(exons) - 1):
                if exons[k][1] == gap_s and exons[k + 1][0] == gap_e:
                    merged_exon = (k, k + 1)
                    break
            new_exon_segments = [(s, e) for k, (s, e, _) in enumerate(exons)]
            if merged_exon is not None:
                k1, k2 = merged_exon
                new_exon_segments = [(s, e) for k, (s, e) in enumerate(new_exon_segments) if k not in (k1, k2)]
                new_exon_segments.append((exons[k1][0], exons[k2][1]))
                new_exon_segments.sort()

            chrom = t['cols'][0]
            strand = t['cols'][6]
            ordered = new_cds_segments if strand == '+' else list(reversed(new_cds_segments))
            seq = cds_sequence(genome, chrom, new_cds_segments, strand)
            protein = translate(seq)
            n_stops = internal_stop_count(protein)

            t['_new_cds'] = new_cds_segments
            t['_new_exon'] = new_exon_segments
            if n_stops == 0:
                fixed_clean.add(tid)
                report.append((tid, 'ShortIntron', f'merged flanking exons ({gap_e - gap_s + 1}bp gap) -> clean ORF'))
            else:
                fixed_still_bad.add(tid)
                report.append((tid, 'ShortIntron', f'merged flanking exons but {n_stops} internal stop(s) remain'))

    # ── Recompute "clean sibling available" now that some ShortIntron fixes succeeded ──
    still_bad = (bad_stop | set(short_intron)) - fixed_clean

    def has_clean_sibling(gid, tid):
        for sib in gene_mrnas.get(gid, []):
            if sib != tid and sib not in still_bad:
                return True
        return False

    drop_transcripts = set()
    pseudogene_transcripts = set()

    for g in genes:
        gid = g['attrs'].get('ID', '')
        for t in g['transcripts']:
            tid = t['old_id']
            if tid not in still_bad:
                continue
            if has_clean_sibling(gid, tid):
                drop_transcripts.add(tid)
                cause = 'StopInProtein' if tid in bad_stop else 'ShortIntron (merge still broken)'
                report.append((tid, cause, 'dropped — clean sibling isoform exists'))
            else:
                pseudogene_transcripts.add(tid)
                cause = 'StopInProtein' if tid in bad_stop else 'ShortIntron (merge still broken)'
                report.append((tid, cause, 'marked pseudogene — sole isoform, no clean alternative'))

    # ── Write output ─────────────────────────────────────────────────────────
    def write_feature(out, cols, a):
        c = cols[:]
        c[8] = attr_str(a) + ';'
        out.write('\t'.join(c) + '\n')

    n_dropped = n_pseudo = n_merged_clean = 0
    with open(args.out, 'w') as out:
        out.write('##gff-version 3\n')
        for g in genes:
            gid = g['attrs'].get('ID', '')
            kept_transcripts = [t for t in g['transcripts'] if t['old_id'] not in drop_transcripts]
            if not kept_transcripts:
                continue   # whole gene had only a dropped transcript — omit gene entirely

            gene_is_pseudo = any(t['old_id'] in pseudogene_transcripts for t in kept_transcripts)
            g_attrs = dict(g['attrs'])
            if gene_is_pseudo:
                g_attrs['pseudo'] = 'true'
                g_attrs['pseudogene'] = 'unknown'
                n_pseudo += 1
            write_feature(out, g['cols'], g_attrs)

            for t in kept_transcripts:
                tid = t['old_id']
                t_attrs = dict(t['attrs'])
                is_pseudo_here = tid in pseudogene_transcripts
                if is_pseudo_here:
                    t_attrs['pseudo'] = 'true'
                write_feature(out, t['cols'], t_attrs)

                if '_new_cds' in t:
                    n_merged_clean += 1 if tid in fixed_clean else 0
                    # emit merged exon list
                    for s, e in sorted(t['_new_exon']):
                        proto = next(ch for ch in t['children'] if ch['type'] == 'exon')
                        ec = proto['cols'][:]
                        ec[3], ec[4] = str(s), str(e)
                        ea = dict(proto['attrs'])
                        write_feature(out, ec, ea)
                    # emit merged CDS list with recomputed phases
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
                    n_dropped += 1

    with open(args.report, 'w') as rf:
        rf.write('transcript_id\tissue\taction\n')
        for tid, issue, action in report:
            rf.write(f"{tid}\t{issue}\t{action}\n")

    import sys
    sys.stderr.write(f"Dropped transcripts (clean sibling existed): {n_dropped}\n")
    sys.stderr.write(f"Marked pseudogene (sole isoform, unfixable): {n_pseudo}\n")
    sys.stderr.write(f"ShortIntron merges that produced a clean ORF: {len(fixed_clean)}\n")
    sys.stderr.write(f"ShortIntron merges still containing a stop: {len(fixed_still_bad)}\n")
    sys.stderr.write(f"Wrote {args.out}\n")
    sys.stderr.write(f"Wrote report: {args.report}\n")


if __name__ == '__main__':
    main()
