#!/usr/bin/env python3
"""Remove only the broken isoforms of genes flagged by funannotate tbl2gbk.

Usage: drop_bad_isoforms.py <genes.gff3> <genome.fasta> <flagged_gene_ids.txt> <out.gff3> [--all-internal-stops]

A transcript of a flagged gene is broken when its spliced CDS has an internal
stop codon or it has an intron (exon or CDS gap) shorter than 10 nt. Broken
transcripts are dropped; the gene is kept with its remaining transcripts and
its span shrunk to fit them. A flagged gene with no remaining transcript, or
in which no transcript is broken by these tests, is removed whole, as before.

With --all-internal-stops the flag file is ignored and EVERY mRNA is tested, for internal
stop codons only (BUSCO-recovery "compleasm" models can carry them and tbl2gbk misses
them). Only the broken mRNAs go. A gene goes only when no mRNA and no other child (such as a
lncRNA) is left; otherwise its span shrinks to the children that stay.
"""
import sys
from collections import defaultdict

gff_file, fasta_file, flag_file, out_file = sys.argv[1:5]
ALL = '--all-internal-stops' in sys.argv[5:]
MIN_INTRON = 10
STOPS = {'TAA', 'TAG', 'TGA'}
_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

flagged = {l.split()[0] for l in open(flag_file) if l.strip()}
if ALL:
    flagged = set()
    for _l in open(gff_file):
        _c = _l.split('\t')
        if len(_c) >= 9 and _c[2] == 'gene': flagged.add(dict(kv.split('=', 1) for kv in _c[8].strip().rstrip(';').split(';') if '=' in kv).get('ID'))

genome, name, chunks = {}, None, []
for line in open(fasta_file):
    if line.startswith('>'):
        if name: genome[name] = ''.join(chunks).upper()
        name, chunks = line[1:].split()[0], []
    else:
        chunks.append(line.strip())
if name: genome[name] = ''.join(chunks).upper()

def attrs(col9):
    return dict(kv.split('=', 1) for kv in col9.strip().rstrip(';').split(';') if '=' in kv)

lines = open(gff_file).read().splitlines()
tx_gene = {}                                  # mRNA id -> gene id (flagged genes only)
exons, cdss = defaultdict(list), defaultdict(list)
for ln in lines:
    c = ln.split('\t')
    if len(c) < 9: continue
    a = attrs(c[8])
    if c[2] in ('mRNA', 'transcript') and a.get('Parent') in flagged:
        tx_gene[a['ID']] = a['Parent']
    elif a.get('Parent') in tx_gene:
        if c[2] == 'exon': exons[a['Parent']].append((int(c[3]), int(c[4])))
        elif c[2] == 'CDS': cdss[a['Parent']].append((c[0], c[6], int(c[3]), int(c[4]), c[7]))

def short_intron(ivs):
    ivs = sorted(ivs)
    return any(b[0] - a[1] - 1 < MIN_INTRON for a, b in zip(ivs, ivs[1:]))

def internal_stop(cds):
    if not cds: return False
    seqid, strand = cds[0][0], cds[0][1]
    parts = sorted(cds, key=lambda x: x[2], reverse=(strand == '-'))
    s = ''.join(genome[seqid][p[2] - 1:p[3]] for p in sorted(cds, key=lambda x: x[2]))
    if strand == '-': s = s.translate(_COMP)[::-1]
    s = s[int(parts[0][4]) if parts[0][4] in '012' else 0:]
    codons = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
    return any(cd in STOPS for cd in codons[:-1])

if ALL:
    bad = {t for t in tx_gene if internal_stop(cdss[t])}
else:
    bad = {t for t in tx_gene
           if short_intron(exons[t]) or short_intron([(x[2], x[3]) for x in cdss[t]]) or internal_stop(cdss[t])}

kept_tx = defaultdict(list)
for t, g in tx_gene.items():
    if t not in bad: kept_tx[g].append(t)
genes_with_bad = {tx_gene[t] for t in bad}
if ALL:
    other_child = set()                       # genes that still have a non-mRNA child (e.g. lncRNA)
    for ln in lines:
        c = ln.split('\t')
        if len(c) >= 9 and c[2] not in ('gene', 'mRNA', 'transcript', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            other_child.add(attrs(c[8]).get('Parent'))
    drop_genes = {g for g in genes_with_bad if not kept_tx[g] and g not in other_child}
else:
    drop_genes = {g for g in flagged if g not in genes_with_bad or not kept_tx[g]}
drop_tx = bad | {t for t, g in tx_gene.items() if g in drop_genes}

span = {}                                     # gene id -> new (start, end) from kept mRNAs
for ln in lines:
    c = ln.split('\t')
    if len(c) >= 9 and c[2] in ('mRNA', 'transcript'):
        a = attrs(c[8])
        if a['ID'] in tx_gene and a['ID'] not in drop_tx:
            s, e = span.get(a['Parent'], (int(c[3]), int(c[4])))
            span[a['Parent']] = (min(s, int(c[3])), max(e, int(c[4])))
if ALL:                                       # re-span only genes that lost an mRNA, over all kept children
    span = {}
    for ln in lines:
        c = ln.split('\t')
        if len(c) >= 9 and c[2] not in ('gene', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            a = attrs(c[8])
            if a.get('Parent') in genes_with_bad and a.get('ID') not in drop_tx:
                s, e = span.get(a['Parent'], (int(c[3]), int(c[4])))
                span[a['Parent']] = (min(s, int(c[3])), max(e, int(c[4])))

with open(out_file, 'w') as out:
    for ln in lines:
        c = ln.split('\t')
        if len(c) >= 9:
            a = attrs(c[8])
            if c[2] == 'gene' and a.get('ID') in drop_genes: continue
            if a.get('ID') in drop_tx or a.get('Parent') in drop_tx: continue
            if c[2] == 'gene' and a.get('ID') in span:
                c[3], c[4] = map(str, span[a['ID']])
                ln = '\t'.join(c)
        out.write(ln + '\n')

print(f'flagged genes: {len(flagged)}; broken isoforms dropped: {len(bad)}; '
      f'genes rescued: {len(flagged) - len(drop_genes)}; genes removed: {len(drop_genes)}'
      if not ALL else f'mRNAs with internal stops dropped: {len(bad)}; genes removed: {len(drop_genes)}', file=sys.stderr)
