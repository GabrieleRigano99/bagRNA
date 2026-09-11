#!/usr/bin/env python3
"""Add non-picked TD2-predicted isoforms to Mikado gene models.

For each complete ORF in the TD2 genome-space GFF3 that is not already in the
filtered pick GFF3, and that:
  1. Overlaps a Mikado-picked gene on the same strand,
  2. Has at least one intron junction confirmed by the portcullis junction BED, and
  3. Has Salmon TPM >= min_tpm for the originating transcript (default 1.0)

…add it as an additional isoform of that gene.

Sources covered: all transcripts in mikado_prepared (Aletsch, StringTie, Trinity,
Helixer, Annevo, Miniprot) for which TD2 predicted a complete ORF.

Portcullis BED12 key: (chrom, thickStart, thickEnd) where
  thickStart = exon1_end (1-based GFF)
  thickEnd   = exon2_start - 1 (1-based GFF)

Chimera guard:
  When an isoform's UTR extends into the genomic territory of a neighboring
  same-strand gene, the overlapping portion is clipped so that the isoform
  stops --boundary-buffer nucleotides before the neighbor starts/ends.
  If the CDS itself would need clipping (true chimeric ORF), the isoform is
  rejected entirely.

Usage:
    add_isoforms.py pick.gff3 td2_genome.gff3 portcullis.bed[,...] output.gff3
                    [--quant-sf quant.sf] [--min-tpm 1.0] [--boundary-buffer 3]
"""
import re
import sys
from collections import defaultdict


# ── helpers ───────────────────────────────────────────────────────────────────

def _strip_orf_suffix(mRNA_id):
    """Strip .pN suffix from TD2 mRNA ID to recover original transcript ID."""
    return re.sub(r'\.p\d+$', '', mRNA_id)


def attr_gff(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


def _splice_sites(exons):
    sites = set()
    for i, (s, e) in enumerate(exons):
        if i > 0:
            sites.add(s)
        if i < len(exons) - 1:
            sites.add(e)
    return sites


def _junctions(exons, chrom):
    result = []
    for i in range(len(exons) - 1):
        result.append((chrom, exons[i][1], exons[i+1][0] - 1))
    return result


# ── loaders ───────────────────────────────────────────────────────────────────

def load_tpm(path):
    """Parse Salmon quant.sf → dict tid: TPM (col1=Name, col4=TPM)."""
    tpm = {}
    with open(path) as f:
        for line in f:
            if line.startswith('Name') or line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 4:
                continue
            try:
                tpm[cols[0]] = float(cols[3])
            except ValueError:
                continue
    return tpm


def load_portcullis(paths):
    """Return frozenset of (chrom, thickStart, thickEnd) from BED12 files."""
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


# ── pick GFF3 parser ──────────────────────────────────────────────────────────

def parse_pick_gff(path):
    """
    Returns:
        genes     : dict gene_id → {chr, strand, start, end, isoforms: list}
        picked_ids: set of original transcript IDs (from alias= attribute)
    """
    genes    = {}
    tx_info  = {}

    with open(path) as f:
        raw = [l for l in f if not l.startswith('#') and l.strip()]

    for line in raw:
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9:
            continue
        feat = cols[2]
        a    = attr_gff(cols[8])
        if feat == 'gene':
            gid = a.get('ID', '')
            genes[gid] = {
                'chr': cols[0], 'strand': cols[6],
                'start': int(cols[3]), 'end': int(cols[4]),
                'isoforms': [],
            }
        elif feat == 'mRNA':
            tid   = a.get('ID', '')
            gid   = a.get('Parent', '')
            alias = a.get('alias', tid)
            tx_info[tid] = {
                'gene_id': gid, 'alias': alias,
                'chr': cols[0], 'strand': cols[6],
                'start': int(cols[3]), 'end': int(cols[4]),
                'exons': [],
            }
        elif feat == 'exon':
            parent = a.get('Parent', '')
            if parent in tx_info:
                tx_info[parent]['exons'].append((int(cols[3]), int(cols[4])))

    picked_ids = set()
    for tid, info in tx_info.items():
        exons = sorted(info['exons'])
        iso = {
            'id': tid, 'chr': info['chr'], 'strand': info['strand'],
            'start': info['start'], 'end': info['end'],
            'splice_sites': _splice_sites(exons),
        }
        gid = info['gene_id']
        if gid in genes:
            genes[gid]['isoforms'].append(iso)
        picked_ids.add(info['alias'])

    return genes, picked_ids


# ── TD2 genome GFF3 parser ────────────────────────────────────────────────────

def parse_td2_genome_gff(path, exclude_ids):
    """
    Parse the TD2 genome-space GFF3 for isoform candidates.

    mRNA IDs are {orig_tid}.pN; orig_tid (without suffix) is used for
    exclusion and TPM lookup. Only multi-exonic ORFs are returned.
    """
    UTR5 = {'five_prime_UTR', '5UTR', 'five_prime_utr'}
    UTR3 = {'three_prime_UTR', '3UTR', 'three_prime_utr'}

    tx       = {}    # td2_mRNA_id → info dict
    skip_set = set() # td2 mRNA IDs whose orig_tid is in exclude_ids

    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feat = cols[2]
            if feat == 'gene':
                continue
            a = attr_gff(cols[8])

            if feat == 'mRNA':
                tid      = a.get('ID', '')
                orig_tid = _strip_orf_suffix(tid)
                if orig_tid in exclude_ids:
                    skip_set.add(tid)
                    continue
                tx[tid] = {
                    'chr': cols[0], 'strand': cols[6],
                    'start': int(cols[3]), 'end': int(cols[4]),
                    'orig_tid': orig_tid, 'source': 'transdecoder',
                    'exons': [], 'cds': [], 'utr5': [], 'utr3': [],
                }

            elif feat == 'exon':
                parent = a.get('Parent', '')
                if parent not in skip_set and parent in tx:
                    tx[parent]['exons'].append((int(cols[3]), int(cols[4])))

            elif feat == 'CDS':
                parent = a.get('Parent', '')
                if parent not in skip_set and parent in tx:
                    phase = int(cols[7]) if cols[7].isdigit() else 0
                    tx[parent]['cds'].append((int(cols[3]), int(cols[4]), phase))

            elif feat in UTR5:
                parent = a.get('Parent', '')
                if parent not in skip_set and parent in tx:
                    tx[parent]['utr5'].append((int(cols[3]), int(cols[4])))

            elif feat in UTR3:
                parent = a.get('Parent', '')
                if parent not in skip_set and parent in tx:
                    tx[parent]['utr3'].append((int(cols[3]), int(cols[4])))

    result = {}
    for tid, info in tx.items():
        if not info['cds']:
            continue
        exons = sorted(info['exons'])
        if len(exons) < 2:
            continue
        info['exons']         = exons
        info['cds']           = sorted(info['cds'])
        info['utr5']          = sorted(info['utr5'])
        info['utr3']          = sorted(info['utr3'])
        info['splice_sites']  = _splice_sites(exons)
        info['junction_keys'] = _junctions(exons, info['chr'])
        result[tid] = info

    return result


# ── gene interval index ───────────────────────────────────────────────────────

def build_gene_interval_index(genes):
    idx = defaultdict(list)
    for gid, g in genes.items():
        idx[g['chr']].append((g['start'], g['end'], g['strand'], gid))
    for chrom in idx:
        idx[chrom].sort()
    return idx


def find_overlapping_gene(idx, chrom, start, end, strand):
    best_gid, best_ov = None, 0
    for (gs, ge, gstrand, gid) in idx.get(chrom, []):
        if ge < start:
            continue
        if gs > end:
            break
        if gstrand != strand:
            continue
        ov = min(end, ge) - max(start, gs) + 1
        if ov > best_ov:
            best_ov, best_gid = ov, gid
    return best_gid


def find_neighbors(idx, chrom, strand, gene_start, gene_end):
    """Return (left_neighbor_end, right_neighbor_start) for same-strand genes flanking
    [gene_start, gene_end]. None where no neighbor exists."""
    left_end = right_start = None
    for (gs, ge, gstrand, gid) in idx.get(chrom, []):
        if gstrand != strand:
            continue
        if ge < gene_start:
            if left_end is None or ge > left_end:
                left_end = ge
        elif gs > gene_end:
            if right_start is None or gs < right_start:
                right_start = gs
    return left_end, right_start


def clip_features(info, clip_start, clip_end):
    """Clip all genomic features of an isoform to [clip_start, clip_end].

    Returns updated info dict, or None if the CDS itself was clipped
    (chimeric ORF) or if fewer than 2 exons remain after clipping.
    """
    # CDS must be fully contained — any clipping means a chimeric ORF
    for (s, e, _) in info['cds']:
        if s < clip_start or e > clip_end:
            return None

    def clip_list(ivals):
        result = []
        for s, e in ivals:
            ns, ne = max(s, clip_start), min(e, clip_end)
            if ns <= ne:
                result.append((ns, ne))
        return result

    new_exons = clip_list(info['exons'])
    if len(new_exons) < 2:
        return None  # became monoexonic — junction check meaningless

    new_info = dict(info)
    new_info['exons'] = sorted(new_exons)
    new_info['utr5']  = clip_list(info['utr5'])
    new_info['utr3']  = clip_list(info['utr3'])
    new_info['start'] = new_exons[0][0]
    new_info['end']   = new_exons[-1][1]
    new_info['splice_sites']  = _splice_sites(new_info['exons'])
    new_info['junction_keys'] = _junctions(new_info['exons'], info['chr'])
    return new_info


# ── GFF3 output ───────────────────────────────────────────────────────────────

def isoform_number(gene):
    return len(gene['isoforms']) + 1


def gff3_lines(gene_id, iso_num, alias, info):
    new_id = f"{gene_id}.{iso_num}"
    chrom, strand, src = info['chr'], info['strand'], info['source']
    lines = [
        '\t'.join([chrom, src, 'mRNA',
                   str(info['start']), str(info['end']), '.', strand, '.',
                   f"ID={new_id};Parent={gene_id};alias={alias};"
                   f"has_start_codon=True;has_stop_codon=True;primary=False;ccode=ab_initio"])
    ]
    for i, (s, e) in enumerate(info['utr5'], 1):
        lines.append('\t'.join([chrom, src, 'five_prime_UTR',
                                str(s), str(e), '.', strand, '.',
                                f"ID={new_id}.utr5p{i};Parent={new_id}"]))
    for i, (s, e) in enumerate(info['exons'], 1):
        lines.append('\t'.join([chrom, src, 'exon',
                                str(s), str(e), '.', strand, '.',
                                f"ID={new_id}.exon{i};Parent={new_id}"]))
    for i, (s, e, ph) in enumerate(info['cds'], 1):
        lines.append('\t'.join([chrom, src, 'CDS',
                                str(s), str(e), '.', strand, str(ph),
                                f"ID=cds.{new_id};Parent={new_id}"]))
    for i, (s, e) in enumerate(info['utr3'], 1):
        lines.append('\t'.join([chrom, src, 'three_prime_UTR',
                                str(s), str(e), '.', strand, '.',
                                f"ID={new_id}.utr3p{i};Parent={new_id}"]))
    return lines


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    raw            = sys.argv[1:]
    quant_sf       = None
    min_tpm        = 1.0
    boundary_buf   = 3
    positional     = []
    i = 0
    while i < len(raw):
        if raw[i] == '--quant-sf' and i + 1 < len(raw):
            quant_sf = raw[i + 1]; i += 2
        elif raw[i] == '--min-tpm' and i + 1 < len(raw):
            min_tpm = float(raw[i + 1]); i += 2
        elif raw[i] == '--boundary-buffer' and i + 1 < len(raw):
            boundary_buf = int(raw[i + 1]); i += 2
        else:
            positional.append(raw[i]); i += 1

    if len(positional) < 4:
        sys.exit(
            f"Usage: {sys.argv[0]} pick.gff3 td2_genome.gff3 portcullis.bed[,...] output.gff3"
            " [--quant-sf quant.sf] [--min-tpm 1.0] [--boundary-buffer 3]"
        )

    pick_path      = positional[0]
    td2_genome_gff = positional[1]
    out_path       = positional[-1]
    port_paths     = positional[2:-1]

    tpm_values = {}
    if quant_sf and quant_sf != 'NO_FILE':
        tpm_values = load_tpm(quant_sf)
        sys.stderr.write(f"Loaded TPM for {len(tpm_values)} transcripts (min_tpm={min_tpm})\n")

    sys.stderr.write("Loading portcullis junctions...\n")
    portcullis = load_portcullis(port_paths)
    sys.stderr.write(f"  {len(portcullis)} verified junctions\n")

    sys.stderr.write("Loading pick GFF3...\n")
    genes, picked_ids = parse_pick_gff(pick_path)
    sys.stderr.write(f"  {len(genes)} genes, {len(picked_ids)} transcripts\n")

    gene_idx = build_gene_interval_index(genes)

    sys.stderr.write("Loading TD2 genome-space ORF models...\n")
    candidates = parse_td2_genome_gff(td2_genome_gff, exclude_ids=picked_ids)
    sys.stderr.write(f"  {len(candidates)} multi-exonic ORFs not in pick\n")

    added       = 0
    added_lines = defaultdict(list)

    for td2_tid, info in sorted(candidates.items()):
        orig_tid = info['orig_tid']

        # 1. Must overlap a Mikado-picked gene on the same strand
        gid = find_overlapping_gene(
            gene_idx, info['chr'], info['start'], info['end'], info['strand']
        )
        if gid is None:
            continue

        # 2. Chimera guard: clip isoform UTRs at midpoint of each intergenic gap.
        #    Using the midpoint (not the neighbor boundary) guarantees that two
        #    isoforms added to adjacent genes cannot overlap each other's UTRs.
        gene_obj = genes[gid]
        left_end, right_start = find_neighbors(
            gene_idx, info['chr'], info['strand'],
            gene_obj['start'], gene_obj['end']
        )
        if left_end is not None:
            clip_s = (left_end + gene_obj['start']) // 2 + boundary_buf
        else:
            clip_s = 1
        if right_start is not None:
            clip_e = (gene_obj['end'] + right_start) // 2 - boundary_buf
        else:
            clip_e = info['end']
        if info['start'] < clip_s or info['end'] > clip_e:
            info = clip_features(info, clip_s, clip_e)
            if info is None:
                continue

        # 3. At least one junction must be portcullis-verified
        if not any(jk in portcullis for jk in info['junction_keys']):
            continue

        # 4. TPM filter — use orig_tid (without .pN) for quant.sf lookup
        if tpm_values and tpm_values.get(orig_tid, 0.0) < min_tpm:
            continue

        gene    = genes[gid]
        iso_num = isoform_number(gene)
        added_lines[gid].extend(gff3_lines(gid, iso_num, orig_tid, info))
        gene['isoforms'].append({
            'id': f"{gid}.{iso_num}",
            'chr': info['chr'], 'strand': info['strand'],
            'start': info['start'], 'end': info['end'],
            'splice_sites': info['splice_sites'],
        })
        added += 1

    sys.stderr.write(f"  Added {added} isoforms to {len(added_lines)} genes\n")

    with open(pick_path) as fin, open(out_path, 'w') as fout:
        current_gene = None
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                fout.write(line)
                continue
            if cols[2] == 'gene':
                if current_gene and current_gene in added_lines:
                    for l in added_lines[current_gene]:
                        fout.write(l + '\n')
                current_gene = attr_gff(cols[8]).get('ID', '')
            fout.write(line)
        if current_gene and current_gene in added_lines:
            for l in added_lines[current_gene]:
                fout.write(l + '\n')


if __name__ == '__main__':
    main()
