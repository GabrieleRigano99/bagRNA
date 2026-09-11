#!/usr/bin/env python3
"""Convert GFF3 + bagRNA annotation TSV to NCBI feature table (.tbl) for table2asn.

Usage: gff3_to_tbl.py <genome.gff3> <functional_annotation.tsv> <out.tbl>

The feature table format spec:
    https://www.ncbi.nlm.nih.gov/genbank/feature_table/

Coordinates are 1-based inclusive.  Minus-strand features are written with
stop < start (first coord > second coord).  Multi-exon features list each
interval on its own line before the qualifiers.
"""
import sys, re
from collections import defaultdict, OrderedDict

gff_file = sys.argv[1]
tsv_file = sys.argv[2]
tbl_file = sys.argv[3]

# ── 1. Parse annotation TSV ───────────────────────────────────────────────
ann = {}   # gene_id → dict of curated fields
with open(tsv_file) as fh:
    header = fh.readline().rstrip('\n').split('\t')
    col = {h: i for i, h in enumerate(header)}
    for line in fh:
        line = line.rstrip('\n')
        if not line:
            continue
        c = line.split('\t')
        if len(c) <= max(col.get('product', 0), col.get('gene_symbol', 0)):
            continue
        gid = c[0]
        ann[gid] = {
            'product':    c[col['product']]         if 'product'           in col else '',
            'symbol':     c[col['gene_symbol']]     if 'gene_symbol'       in col else '',
            'ec':         [e for e in c[col['EC_numbers']].split(',')  if e] if 'EC_numbers'       in col else [],
            'go':         [g for g in c[col['GO_terms']].split('|')   if g] if 'GO_terms'          in col else [],
            'interpro':   [i for i in c[col['InterPro_accessions']].split('|') if i] if 'InterPro_accessions' in col else [],
            'pfam':       [p for p in c[col['Pfam_domains']].split('|') if p] if 'Pfam_domains'   in col else [],
            'kegg_ko':    [k for k in c[col['KEGG_KO']].split(',')    if k] if 'KEGG_KO'          in col else [],
            'cog':        c[col['COG_category']]    if 'COG_category'      in col else '',
            'nog_og':     c[col['eggNOG_OG']]       if 'eggNOG_OG'        in col else '',
            'signalp':    c[col['SignalP']]  == 'Y' if 'SignalP'           in col else False,
            'secreted':   c[col['Secreted']] == 'Y' if 'Secreted'          in col else False,
            'smcog':      c[col['SMCOG']]           if 'SMCOG'             in col else '',
            'bgc_type':   c[col['BGC_cluster_type']] if 'BGC_cluster_type' in col else '',
            'bgc_role':   c[col['BGC_gene_role']]   if 'BGC_gene_role'    in col else '',
            'effector':   c[col['EffectorP_class']] if 'EffectorP_class'  in col else '',
            'cazyme':     c[col['CAZyme_family']]   if 'CAZyme_family'    in col else '',
            'merops':     c[col['MEROPS_family']]   if 'MEROPS_family'    in col else '',
            'phi_acc':    c[col['PHIbase_accession']] if 'PHIbase_accession' in col else '',
        }

def best_product(gid):
    a = ann.get(gid, {})
    p = a.get('product', '').strip()
    if not p or p.lower() in ('hypothetical protein', ''):
        return 'hypothetical protein'
    # Strip leading qualifiers that NCBI rejects
    p = re.sub(r'^(putative|probable|possible|predicted|uncharacterized|unknown)\s+', '', p, flags=re.I)
    p = p.rstrip('.')
    return p if p else 'hypothetical protein'

# ── 2. Parse GFF3 ────────────────────────────────────────────────────────
# Build per-sequence ordered gene list with full subfeature coordinates.

class Gene:
    __slots__ = ('seqid', 'start', 'stop', 'strand', 'attrs',
                 'locus_tag', 'symbol', 'mrnas')
    def __init__(self):
        self.mrnas = OrderedDict()   # mrna_id → MRna

class MRna:
    __slots__ = ('start', 'stop', 'strand', 'exons', 'cdss')
    def __init__(self):
        self.exons = []   # (start, stop)
        self.cdss  = []   # (start, stop)

genes    = OrderedDict()     # gene_id → Gene
mrna_map = {}                # mrna_id → gene_id
seq_genes = defaultdict(list)  # seqid → [gene_id, ...]

def parse_attrs(attr_str):
    d = {}
    for kv in attr_str.rstrip(';').split(';'):
        if '=' in kv:
            k, _, v = kv.partition('=')
            d[k.strip()] = v.strip()
    return d

with open(gff_file) as fh:
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9:
            continue
        seqid, _, ftype, start, stop, _, strand, _, attr_str = cols[:9]
        start, stop = int(start), int(stop)
        attrs = parse_attrs(attr_str)
        gid  = attrs.get('ID', '')
        par  = attrs.get('Parent', '')

        if ftype == 'gene' and gid:
            g = Gene()
            g.seqid    = seqid
            g.start    = start
            g.stop     = stop
            g.strand   = strand
            g.locus_tag = attrs.get('locus_tag', gid)
            g.symbol   = attrs.get('gene', attrs.get('Name', ''))
            genes[gid] = g
            seq_genes[seqid].append(gid)

        elif ftype in ('mRNA', 'transcript') and gid and par:
            if par not in genes:
                continue
            m = MRna()
            m.start  = start
            m.stop   = stop
            m.strand = strand
            genes[par].mrnas[gid] = m
            mrna_map[gid] = par

        elif ftype == 'exon' and par:
            gene_id = mrna_map.get(par)
            if gene_id and par in genes[gene_id].mrnas:
                genes[gene_id].mrnas[par].exons.append((start, stop))

        elif ftype == 'CDS' and par:
            gene_id = mrna_map.get(par)
            if gene_id and par in genes[gene_id].mrnas:
                genes[gene_id].mrnas[par].cdss.append((start, stop))

# ── 3. Write .tbl ────────────────────────────────────────────────────────

def tbl_coords(intervals, strand):
    """Return list of (c1, c2) in feature-table orientation."""
    ivs = sorted(intervals, key=lambda x: x[0])
    if strand == '-':
        return [(s2, s1) for s1, s2 in ivs]
    return ivs

def write_feature(out, ftype, intervals, strand, qualifiers):
    ivs = tbl_coords(intervals, strand)
    for i, (c1, c2) in enumerate(ivs):
        if i == 0:
            out.write(f'{c1}\t{c2}\t{ftype}\n')
        else:
            out.write(f'{c1}\t{c2}\n')
    for key, val in qualifiers:
        if val:
            out.write(f'\t\t\t{key}\t{val}\n')

def build_note(gid):
    a = ann.get(gid, {})
    parts = []
    if a.get('smcog'):   parts.append('smCOG: ' + a['smcog'])
    if a.get('bgc_type'): parts.append('Secondary metabolite cluster: ' + a['bgc_type'])
    elif a.get('bgc_role') and a['bgc_role'] != 'other':
        parts.append('BGC gene (' + a['bgc_role'] + ')')
    if a.get('effector'):  parts.append('Predicted ' + a['effector'].lower() + ' effector (EffectorP3)')
    if a.get('cazyme'):    parts.append('CAZyme: ' + a['cazyme'])
    if a.get('merops'):    parts.append('Peptidase family ' + a['merops'] + ' (MEROPS)')
    if a.get('secreted'):  parts.append('Predicted secreted protein')
    return '; '.join(parts)

with open(tbl_file, 'w') as out:
    for seqid, gid_list in seq_genes.items():
        out.write(f'>Feature {seqid}\n')

        for gid in gid_list:
            g = genes[gid]
            a = ann.get(gid, {})
            product  = best_product(gid)
            symbol   = a.get('symbol', '') or g.symbol or ''
            locus    = g.locus_tag

            # gene feature — single span
            gene_quals = [('locus_tag', locus)]
            if symbol:
                gene_quals.append(('gene', symbol))
            write_feature(out, 'gene', [(g.start, g.stop)], g.strand, gene_quals)

            for mid, m in g.mrnas.items():
                exon_ivs = m.exons if m.exons else [(m.start, m.stop)]
                cds_ivs  = m.cdss  if m.cdss  else []

                # mRNA feature
                mrna_quals = [('locus_tag', locus), ('product', product)]
                if symbol:
                    mrna_quals.insert(0, ('gene', symbol))
                write_feature(out, 'mRNA', exon_ivs, m.strand, mrna_quals)

                # CDS feature
                if cds_ivs:
                    cds_quals = []
                    if symbol:
                        cds_quals.append(('gene', symbol))
                    cds_quals.append(('locus_tag', locus))
                    cds_quals.append(('product', product))
                    cds_quals.append(('codon_start', '1'))
                    for ec in a.get('ec', []):
                        cds_quals.append(('EC_number', ec))
                    for go in a.get('go', []):
                        cds_quals.append(('db_xref', 'GO:' + go if not go.startswith('GO:') else go))
                    for ip in a.get('interpro', []):
                        cds_quals.append(('db_xref', 'InterPro:' + ip if not ip.startswith('IPR') else 'InterPro:' + ip))
                    for ko in a.get('kegg_ko', []):
                        cds_quals.append(('db_xref', 'KEGG:' + ko if not ko.startswith('K') else 'KEGG:' + ko))
                    note = build_note(gid)
                    if note:
                        cds_quals.append(('note', note))
                    write_feature(out, 'CDS', cds_ivs, m.strand, cds_quals)

print(f'Wrote {tbl_file}: {sum(len(v) for v in seq_genes.values())} genes across {len(seq_genes)} sequences')
