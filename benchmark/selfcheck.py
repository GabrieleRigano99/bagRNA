#!/usr/bin/env python3
"""Tiny self-check for bin/gff3_to_tbl.py and bin/drop_bad_isoforms.py (stdlib only). Run: python3 benchmark/selfcheck.py"""
import os, subprocess, sys, tempfile
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin')
d = tempfile.mkdtemp()
def w(name, text):
    open(os.path.join(d, name), 'w').write(text); return os.path.join(d, name)
def run(script, *a): return subprocess.run([sys.executable, '-I', os.path.join(BIN, script), *a], capture_output=True, text=True, check=True)

# contig: 40 bp ATG..TAA gene (+), 20 N gap, 2nd gene (+) ends in the gap, 3rd gene with an internal stop
orf = 'ATG' + 'GCT' * 10 + 'TAA'
seq = orf + 'C' * 10 + 'N' * 20 + 'C' * 10 + orf + 'C' * 10 + 'ATG' + 'TAA' + 'GCT' * 3 + 'TAG'
fa = w('g.fa', '>c1\n' + seq + '\n')
L = len(orf); s2 = L + 10 + 20 + 10 + 1; s3 = s2 + L + 10
def gene(i, a, b, extra=''):
    return (f'c1\tx\tgene\t{a}\t{b}\t.\t+\t.\tID=G{i};locus_tag=T_{i}{extra}\n'
            f'c1\tx\tmRNA\t{a}\t{b}\t.\t+\t.\tID=G{i}-T1;Parent=G{i}\n'
            f'c1\tx\texon\t{a}\t{b}\t.\t+\t.\tID=G{i}-T1.e;Parent=G{i}-T1\n'
            f'c1\tx\tCDS\t{a}\t{b}\t.\t+\t0\tID=G{i}-T1.c;Parent=G{i}-T1\n')
# gene 1 and 2 share symbol "ABC"; gene 4 ends inside the N gap
gff = w('a.gff3', gene(1, 1, L, ';gene=ABC') + gene(2, s2, s2 + L - 1, ';gene=ABC') + gene(3, s3, s3 + 14) + gene(4, L - 5, L + 15) + gene(5, 2, L))
tsv = w('a.tsv', 'gene_id\tproduct\tgene_symbol\n' + ''.join(f'G{i}\tkinase\t\n' for i in (1, 2, 3, 4, 5)))
out = os.path.join(d, 'o.tbl'); run('gff3_to_tbl.py', gff, tsv, out, fa)
tbl = open(out).read()
assert '\t\t\tgene\tABC' not in tbl, 'shared symbol must not be written as /gene'
assert 'T_4' not in tbl, 'gene ending in an N gap must be left out'
assert 'T_1' in tbl and 'T_2' in tbl and 'T_3' in tbl
blk5 = tbl.split('\tgene\n')[[i for i, b in enumerate(tbl.split('\tgene\n')) if 'T_5' in b][0]]
assert '>' in blk5, 'CDS of length not a multiple of 3 must be marked 3-prime partial (NoStop)'

keep = os.path.join(d, 'k.gff3'); open(os.path.join(d, 'none.txt'), 'w').close()
run('drop_bad_isoforms.py', gff, fa, os.path.join(d, 'none.txt'), keep, '--all-internal-stops')
k = open(keep).read()
assert 'G3' not in k, 'internal-stop mRNA (and its childless gene) must be dropped'
assert 'ID=G1-T1' in k and 'ID=G2-T1' in k, 'healthy mRNAs must stay'
print('selfcheck OK')
