#!/usr/bin/env python3
"""Reformat agat_sp_manage_IDs.pl --tair --prefix output into NCBI-style
locus-tag IDs, matching the convention this pipeline already relies on
(and that funannotate util gff2tbl/tbl2gbk expect verbatim, unreformatted,
in the GFF's own ID= attribute):

    gene              PREFIX_NNNNNN            (6-digit zero-padded)
    mRNA/tRNA/ncRNA/rRNA (isoform M of gene N)  PREFIX_NNNNNN-TM
    exon (segment K)                            PREFIX_NNNNNN-TM.exonK
    five_prime_UTR (segment K)                  PREFIX_NNNNNN-TM.utr5pK
    three_prime_UTR (segment K)                 PREFIX_NNNNNN-TM.utr3pK
    CDS (all segments share one ID)             PREFIX_NNNNNN-TM.cds

AGAT's --tair mode already assigns clean, hierarchical, deterministic IDs
(PREFIXn / PREFIXn.m / PREFIXn.m-<kind><k>) and — unlike funannotate
gff-rename — never touches any other attribute (product=, isotype=,
anticodon=, Name=, Alias= all pass through untouched). This script
rewrites ID=/Parent= syntax to match (it does not re-derive the
hierarchy), sets a real product= for tRNA (from the isotype= tRNAscan-SE
attribute, still present here — unlike after funannotate gff-rename) and
ncRNA (hypothetical lncRNA/ncRNA by length, matching this pipeline's
existing convention), then strips Alias=/Name= for a clean final GFF.
funannotate util gff2tbl synthesises "product None" on its own for
anything left unset, so this only needs to cover the cases we have real
data for.

Usage:
    reformat_locus_tag_ids.py input.gff3 PREFIX output.gff3
"""
import re
import sys

SUFFIX_MAP = {
    'exon': ('exon', True),
    'cds': ('cds', False),
    'five_prime_utr': ('utr5p', True),
    'three_prime_utr': ('utr3p', True),
}


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} input.gff3 PREFIX output.gff3")

    in_path, prefix, out_path = sys.argv[1:4]
    pfx = re.escape(prefix)

    gene_re = re.compile(rf'^{pfx}(\d+)$')
    level2_re = re.compile(rf'^{pfx}(\d+)\.(\d+)$')
    child_re = re.compile(rf'^{pfx}(\d+)\.(\d+)-([a-zA-Z_]+)(\d+)$')

    def new_gene_id(n):
        return f"{prefix}_{int(n):06d}"

    def new_level2_id(n, m):
        return f"{new_gene_id(n)}-T{m}"

    def remap(token):
        """Rewrite a single ID/Parent value if it matches an AGAT --tair pattern."""
        m = child_re.match(token)
        if m:
            n, iso, kind, k = m.groups()
            base = new_level2_id(n, iso)
            suffix, per_segment_unique = SUFFIX_MAP.get(kind.lower(), (kind, True))
            return f"{base}.{suffix}{k}" if per_segment_unique else f"{base}.{suffix}"
        m = level2_re.match(token)
        if m:
            n, iso = m.groups()
            return new_level2_id(n, iso)
        m = gene_re.match(token)
        if m:
            return new_gene_id(m.group(1))
        return token  # not an AGAT-assigned ID (e.g. unrelated value) — leave as-is

    def attr_dict(field):
        d = {}
        for tok in field.split(';'):
            tok = tok.strip()
            if tok.startswith('ID=') or tok.startswith('Parent='):
                continue  # handled separately, order-sensitive
            if '=' in tok:
                k, v = tok.split('=', 1)
                d[k] = v
        return d

    n_lines = 0
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                fout.write(line)
                continue

            feat = cols[2]
            parts = cols[8].split(';')
            for i, tok in enumerate(parts):
                if tok.startswith('ID='):
                    parts[i] = 'ID=' + remap(tok[3:])
                elif tok.startswith('Parent='):
                    parts[i] = 'Parent=' + remap(tok[7:])

            # Idempotent: this script runs twice (before and after tbl2gbk QC
            # removal), and the second pass sees lines the first pass already
            # tagged — never append a second product=/ncRNA_class= attribute.
            existing = attr_dict(cols[8])
            if feat == 'tRNA' and 'product' not in existing:
                isotype = existing.get('isotype')
                if isotype and isotype.lower() != 'undet':
                    parts.append(f'product=tRNA-{isotype}')
            elif feat == 'ncRNA' and 'product' not in existing:
                length = int(cols[4]) - int(cols[3]) + 1
                if length >= 200:
                    parts.append('product=hypothetical lncRNA')
                    parts.append('ncRNA_class=lncRNA')
                else:
                    parts.append('product=hypothetical ncRNA')
                    parts.append('ncRNA_class=other')

            parts = [p for p in parts if not p.startswith('Alias=') and not p.startswith('Name=')]
            cols[8] = ';'.join(parts)
            fout.write('\t'.join(cols) + '\n')
            n_lines += 1

    sys.stderr.write(f"Reformatted {n_lines} feature lines to locus-tag IDs.\n")


if __name__ == '__main__':
    main()
