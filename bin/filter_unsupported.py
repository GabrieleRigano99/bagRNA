#!/usr/bin/env python3
"""Report multi-exonic ab initio gene models with no portcullis-confirmed junction.

For Helixer/Annevo primary picks where none of the predicted introns are supported
by RNA-seq splice junctions, the model is left untouched here — collapsing it to a
single exon would merge former intron sequence into the CDS (former phases become
meaningless and the frame is corrupted for everything past each removed intron),
which reliably manufactures internal stop codons. FINAL_JUNCTION_FILTER (run later,
after isoforms are added) already makes the correct call per Mikado's own logic: a
multi-exonic model with unconfirmed junctions produces a better protein than a
collapsed monoexonic one, so unsupported models are only dropped there if the gene
still has a supported alternative; otherwise they are kept as-is.

Protein-evidence models (Miniprot / ProteinEvidence) are never flagged since
they are supported by protein homology rather than RNA-seq.
Monoexonic models are never flagged.

Usage:
    filter_unsupported.py pick.gff3 portcullis.bed[,...] output.gff3
"""
import shutil
import sys


def load_portcullis(paths):
    junctions = set()
    for path in paths:
        with open(path) as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                c = line.split('\t')
                if len(c) < 8:
                    continue
                try:
                    junctions.add((c[0], int(c[6]), int(c[7])))
                except (ValueError, IndexError):
                    continue
    return frozenset(junctions)


def attr_gff(field):
    d = {}
    for tok in field.split(';'):
        tok = tok.strip()
        if '=' in tok:
            k, v = tok.split('=', 1)
            d[k] = v
    return d


def infer_source_from_alias(alias, col_source):
    al = alias.lower()
    if 'proteinevidence' in al or 'miniprot' in al:
        return 'ProteinEvidence'
    if 'helixer' in al:
        return 'Helixer'
    if 'annevo' in al:
        return 'Annevo'
    if 'stringtie' in al:
        return 'StringTie'
    if 'aletsch' in al:
        return 'Aletsch'
    if 'trinity' in al:
        return 'Trinity'
    return col_source


def main():
    if len(sys.argv) < 4:
        sys.exit(f"Usage: {sys.argv[0]} pick.gff3 portcullis.bed[,...] output.gff3")

    pick_path  = sys.argv[1]
    out_path   = sys.argv[-1]
    port_paths = sys.argv[2:-1]

    sys.stderr.write("Loading portcullis junctions...\n")
    portcullis = load_portcullis(port_paths)
    sys.stderr.write(f"  {len(portcullis)} verified junctions\n")

    # ── Collect mRNA metadata and exons (for reporting only) ─────────────────
    mrnas   = {}   # tid → metadata dict
    cur_tid = None

    with open(pick_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feat = cols[2]
            a = attr_gff(cols[8])

            if feat == 'mRNA':
                cur_tid = a.get('ID', '')
                alias   = a.get('alias', cur_tid)
                mrnas[cur_tid] = {
                    'chr': cols[0],
                    'source': infer_source_from_alias(alias, cols[1]),
                    'exons': [],
                }
            elif feat == 'exon' and cur_tid:
                parent = a.get('Parent', cur_tid)
                if parent in mrnas:
                    mrnas[parent]['exons'].append((int(cols[3]), int(cols[4])))
            elif feat == 'gene':
                cur_tid = None

    PROTEIN_SOURCES = {'ProteinEvidence', 'Miniprot'}
    n_unsupported = 0

    for info in mrnas.values():
        exons = sorted(info['exons'])
        if len(exons) < 2 or info['source'] in PROTEIN_SOURCES:
            continue
        jkeys = [(info['chr'], exons[i][1], exons[i + 1][0] - 1)
                 for i in range(len(exons) - 1)]
        if not any(j in portcullis for j in jkeys):
            n_unsupported += 1

    sys.stderr.write(
        f"{n_unsupported} ab initio multi-exonic model(s) have no portcullis-"
        f"confirmed junction; left untouched (FINAL_JUNCTION_FILTER decides "
        f"their fate downstream)\n"
    )

    shutil.copyfile(pick_path, out_path)


if __name__ == '__main__':
    main()
