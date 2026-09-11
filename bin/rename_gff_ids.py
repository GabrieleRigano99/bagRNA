#!/usr/bin/env python3
"""Uniformly rename all gene/isoform IDs in the final structural GFF3.

Replicates funannotate gff-rename's convention: genes numbered sequentially
in genomic order as PREFIX_NNNNNN, isoforms as PREFIX_NNNNNN-T1/-T2/-T3...
(preserving each gene's existing child order), children as
{mrna_id}.exonN / {mrna_id}.cds / {mrna_id}.utr5pN / {mrna_id}.utr3pN.

Why this is needed: struct_final_v2.gff3 mixes several ID schemes because
different bin/ scripts added loci at different times —
  - original pipeline genes:      SS02_000001 / SS02_000001-T1
  - TD2-added isoforms:           SS02_000001.2, SS02_000001.3
  - TD2 standalone new loci:      raw TransDecoder IDs leaked through
                                   unchanged, e.g.
                                   "Aletsch_MSTRG.3086.1.gene^CBS145945_2^+"
  - BUSCO-recovery isoforms:      busco_iso_000001 (no gene prefix at all)
  - BUSCO standalone new genes:   busco_gene_r000001
This script produces one uniform scheme across all of them.

Also strips pipeline-internal bookkeeping attributes (alias, ccode,
has_start_codon, has_stop_codon, primary, Name) and the hardcoded
placeholder `product=` on mRNA/ncRNA (always "hypothetical protein" /
"hypothetical lncRNA" pre-rename) so a downstream functional-annotation
merge can actually set real values — merge_attrs()-style mergers skip any
key already present, so a placeholder left in place silently blocks every
real annotation from ever being written.

Usage:
    rename_gff_ids.py input.gff3 output.gff3 id_map.tsv [--prefix SS02]
"""
import sys
import argparse
from collections import defaultdict

STRIP_KEYS = {'alias', 'Alias', 'ccode', 'has_start_codon', 'has_stop_codon', 'primary', 'Name',
              # tRNA anticodon= (from raw tRNAscan-SE, kept on recovered loci): table2asn
              # expects a specific /anticodon=(pos:...,aa:...,seq:...) qualifier it can't
              # derive from a bare GFF3 attribute, and flags it as SEQ_FEAT.UnparsedtRNAAnticodon.
              # The 78 originally-surviving tRNA never carried this attribute either.
              'anticodon'}
TRANSCRIPT_TYPES = {'mRNA', 'ncRNA', 'tRNA', 'rRNA', 'lncRNA', 'transcript'}


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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input_gff')
    p.add_argument('output_gff')
    p.add_argument('id_map_tsv')
    p.add_argument('--prefix', default='SS02')
    p.add_argument('--trnascan-gff',
                    help="Raw tRNAscan-SE --gff output (e.g. run10july/trnascan/trnascan.gff). "
                         "funannotate gff-rename stringifies tRNA features with no /product "
                         "qualifier as the literal text 'product=None' — reconstruct the real "
                         "product ('tRNA-Ala', etc.) from the isotype= attribute by matching "
                         "tRNA coordinates back to this file.")
    return p.parse_args()


def load_trna_isotypes(path):
    """(chrom, start, end, strand) -> isotype, from raw tRNAscan-SE --gff output."""
    lut = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9 or c[2] != 'tRNA':
                continue
            m = attr(c[8]).get('isotype')
            if m:
                lut[(c[0], int(c[3]), int(c[4]), c[6])] = m
    return lut


def main():
    args = parse_args()
    trna_isotypes = load_trna_isotypes(args.trnascan_gff) if args.trnascan_gff else {}

    # ── Parse into gene blocks, preserving each gene's contiguous children ──
    chrom_order = []
    seen_chroms = set()
    genes = []   # list of dicts: {chrom, start, end, strand, score, line, children: [...]}
    cur_gene = None

    with open(args.input_gff) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            feat = c[2]
            a = attr(c[8])

            if feat == 'gene':
                gid = a.get('ID', '')
                if c[0] not in seen_chroms:
                    seen_chroms.add(c[0])
                    chrom_order.append(c[0])
                cur_gene = {
                    'old_id': gid, 'chrom': c[0], 'start': int(c[3]), 'end': int(c[4]),
                    'strand': c[6], 'score': c[5], 'source': c[1],
                    'transcripts': [],   # list of {old_id, type, cols, attrs, children:[...]}
                }
                genes.append(cur_gene)
            elif feat in TRANSCRIPT_TYPES:
                if cur_gene is None:
                    continue
                tid = a.get('ID', '')
                t = {'old_id': tid, 'type': feat, 'cols': c[:], 'attrs': a, 'children': []}
                cur_gene['transcripts'].append(t)
            else:
                if cur_gene is None or not cur_gene['transcripts']:
                    continue
                parent = a.get('Parent', '')
                # attach to whichever transcript in the current gene owns this Parent
                for t in reversed(cur_gene['transcripts']):
                    if t['old_id'] == parent:
                        t['children'].append({'type': feat, 'cols': c[:], 'attrs': a})
                        break

    chrom_rank = {c: i for i, c in enumerate(chrom_order)}
    genes.sort(key=lambda g: (chrom_rank[g['chrom']], g['start']))

    id_map = []   # (old_id, new_id) pairs — genes and transcripts only
    out_lines = ['##gff-version 3\n']

    for gi, gene in enumerate(genes, 1):
        new_gid = f"{args.prefix}_{gi:06d}"
        id_map.append((gene['old_id'], new_gid))

        gene_attrs = {'ID': new_gid}
        out_lines.append('\t'.join([
            gene['chrom'], gene['source'], 'gene',
            str(gene['start']), str(gene['end']), gene['score'], gene['strand'], '.',
            attr_str(gene_attrs) + ';',
        ]) + '\n')

        for ti, t in enumerate(gene['transcripts'], 1):
            new_tid = f"{new_gid}-T{ti}"
            id_map.append((t['old_id'], new_tid))

            ta = {k: v for k, v in t['attrs'].items() if k not in STRIP_KEYS}
            ta['ID'] = new_tid
            ta['Parent'] = new_gid
            if t['type'] in ('ncRNA', 'mRNA') and 'product' in ta:
                del ta['product']   # placeholder — let functional merge set the real value
            if t['type'] == 'tRNA':
                # funannotate gff-rename stringifies tRNA with no /product
                # qualifier as literal "product=None" — tRNAscan-SE's own
                # output never had a product field at all, only isotype=.
                # Always (re)derive from the coordinate-matched isotype when
                # available, whether the existing value is "None", missing
                # entirely (recovered loci), or anything else.
                key = (t['cols'][0], int(t['cols'][3]), int(t['cols'][4]), t['cols'][6])
                isotype = trna_isotypes.get(key)
                if isotype:
                    ta['product'] = f"tRNA-{isotype}"
                elif ta.get('product') == 'None':
                    ta['product'] = 'tRNA'

            tc = t['cols'][:]
            tc[8] = attr_str(ta) + ';'
            out_lines.append('\t'.join(tc) + '\n')

            exon_n = 0
            utr5_n = 0
            utr3_n = 0
            for ch in t['children']:
                ca = {k: v for k, v in ch['attrs'].items() if k not in STRIP_KEYS}
                ca['Parent'] = new_tid
                if ch['type'] == 'exon':
                    exon_n += 1
                    ca['ID'] = f"{new_tid}.exon{exon_n}"
                elif ch['type'] == 'CDS':
                    ca['ID'] = f"{new_tid}.cds"
                elif ch['type'] == 'five_prime_UTR':
                    utr5_n += 1
                    ca['ID'] = f"{new_tid}.utr5p{utr5_n}"
                elif ch['type'] == 'three_prime_UTR':
                    utr3_n += 1
                    ca['ID'] = f"{new_tid}.utr3p{utr3_n}"
                else:
                    ca['ID'] = f"{new_tid}.{ch['type']}"
                cc = ch['cols'][:]
                cc[8] = attr_str(ca) + ';'
                out_lines.append('\t'.join(cc) + '\n')

    with open(args.output_gff, 'w') as out:
        out.writelines(out_lines)

    with open(args.id_map_tsv, 'w') as out:
        out.write('old_id\tnew_id\n')
        for old, new in id_map:
            out.write(f"{old}\t{new}\n")

    sys.stderr.write(f"Renamed {len(genes)} genes, {len(id_map) - len(genes)} transcripts\n")
    sys.stderr.write(f"Wrote {args.output_gff}\n")
    sys.stderr.write(f"Wrote ID map: {args.id_map_tsv}\n")


if __name__ == '__main__':
    main()
