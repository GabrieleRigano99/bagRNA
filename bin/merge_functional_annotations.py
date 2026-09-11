#!/usr/bin/env python3
"""Merge all functional annotation tool outputs into the final structural GFF3.

Standalone replacement for modules/annotate_functional.nf's embedded merger.
Fixes three bugs found in that version:
  - eggNOG: column-count guard (`len(c) < 22`) never matched the real
    21-column emapper output, so every eggNOG row was silently skipped.
  - EffectorP3: exact-matched 'Y' against the real "Y (0.92)" probability
    strings, so no effector call ever matched.
  - TMbed: read the descriptor from the wrong column (interpro TSV col 5,
    always '-') instead of col 4, where "TMhelix_*"/"TMbeta_*" actually is.

Annotations are aggregated at the gene level (shared across all isoforms of
a gene) and written as GFF3 attributes on gene/mRNA/CDS features, plus a
flat summary TSV and a stats report.

Usage:
    merge_functional_annotations.py --gff final.gff3 --indir functional_annotation/ \\
        --out-gff annotated.gff3 --out-tsv functional_annotation.tsv --out-stats annotation_stats.txt
"""
import argparse
import os
import re
from collections import defaultdict


def is_absent(f):
    return not f or not os.path.exists(f) or os.path.getsize(f) == 0


def new_entry():
    return {
        'mrna_ids': set(), 'product': '', 'gene_symbol': '',
        'go': set(), 'go_propagated': set(), 'ec': set(),
        'kegg_ko': set(), 'kegg_path': set(), 'cog': '', 'ko_primary': '',
        'nog_og': '', 'nog_desc': '',
        'interpro': set(), 'pfam': set(), 'panther': '',
        'signalp': False, 'tm_tmbed': 0, 'tm_phobius': 0, 'sp_phobius': False,
        'effector': '', 'cazyme': '', 'caz_tools': 0, 'pul': '',
        'mer_hit': '', 'mer_fam': '', 'tc': '',
        'phi_acc': '', 'phi_gene': '', 'phi_pheno': '',
        'rfam': set(),
        'bgc_type': '', 'bgc_role': '', 'bgc_domains': set(), 'smcog': '',
    }


_STRIP_LEAD = re.compile(r'^(putative|probable|possible|predicted|uncharacterized|unknown)\s+', re.I)
_BAD_PRODUCT = {
    '', 'hypothetical protein', 'uncharacterized protein', 'unknown protein',
    'expressed protein', 'unknown function', 'unnamed protein product',
}


def ncbi_product(name, symbol=''):
    if not name:
        return 'hypothetical protein'
    name = _STRIP_LEAD.sub('', name.strip()).rstrip('.')
    # A product that's just the gene symbol repeated (e.g. product="APC11",
    # gene_symbol="APC11") isn't a description — some upstream source gave
    # us a bare symbol where a description was expected. Fall back rather
    # than display a symbol as if it explained what the gene product is.
    if symbol and name.lower() == symbol.strip().lower():
        return 'hypothetical protein'
    return name if name.lower() not in _BAD_PRODUCT else 'hypothetical protein'


def is_secreted(a):
    return (a['signalp'] or a['sp_phobius']) and a['tm_tmbed'] == 0 and a['tm_phobius'] == 0


def gff_val(s):
    return s.replace('%', '%25').replace(';', '%3B').replace('=', '%3D').replace(',', '%2C')


def merge_attrs(attrs_str, extras):
    # The functional annotation is authoritative for the keys it sets
    # (product, gene, Name, EC_number, Dbxref, Ontology_term, Note). Drop any
    # pre-existing occurrence of those keys from the structural GFF and replace
    # with our value — otherwise a placeholder baked in upstream (e.g. the
    # pipeline's "product=hypothetical protein" on every mRNA) would block the
    # real product. Keys we don't set (ID, Parent, ncRNA_class, ...) are kept
    # in their original order.
    set_keys = {k for k, v in extras if v}
    kept = []
    for kv in attrs_str.rstrip(';').split(';'):
        if not kv:
            continue
        k = kv.split('=', 1)[0] if '=' in kv else kv
        if k not in set_keys:
            kept.append(kv)
    result = ';'.join(kept)
    for key, val in extras:
        if val:
            result = (result + ';' if result else '') + key + '=' + val
    return result


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gff', required=True, help='Final structural annotation GFF3')
    p.add_argument('--indir', help='functional_annotation/ output dir; used to fill in default file paths below')
    p.add_argument('--eggnog')
    p.add_argument('--interpro')
    p.add_argument('--kofamscan')
    p.add_argument('--ko-info', help='KO-keyed ko_info.tsv from bin/fetch_kegg_ko_info.py — KEGG REST symbol/name '
                                      'lookup covering KOs from both KofamScan and eggNOG')
    p.add_argument('--infernal')
    p.add_argument('--effectorp3')
    p.add_argument('--dbcan-overview')
    p.add_argument('--dbcan-cgc')
    p.add_argument('--merops')
    p.add_argument('--phi-base')
    p.add_argument('--dbcan-tc', help='dbCAN diamond.out.tc — TCDB transporter classification hits')
    p.add_argument('--go-obo')
    p.add_argument('--antismash-gbk')
    p.add_argument('--gene2product')
    p.add_argument('--id-map', help='old_id\\tnew_id TSV from rename_gff_ids.py — translates tool-output '
                                     'query IDs (from before renaming) to the current GFF ID namespace')
    p.add_argument('--out-gff', required=True)
    p.add_argument('--out-tsv', required=True)
    p.add_argument('--out-stats', required=True)
    args = p.parse_args()

    if args.indir:
        def d(explicit, relpath):
            return explicit or os.path.join(args.indir, relpath)
        args.eggnog         = d(args.eggnog,         'eggnog_output.emapper.annotations')
        args.interpro       = d(args.interpro,       'interpro_output.tsv')
        args.kofamscan      = d(args.kofamscan,      'kofamscan_result.tsv')
        args.ko_info = d(args.ko_info, 'ko_info.tsv')
        args.infernal       = d(args.infernal,       'infernal_table.txt')
        args.effectorp3     = d(args.effectorp3,     'effectorp3_output.txt')
        args.dbcan_overview = d(args.dbcan_overview, 'dbcan/overview.tsv')
        args.dbcan_cgc      = d(args.dbcan_cgc,      'dbcan/cgc_standard_out.tsv')
        args.merops         = d(args.merops,         'merops/merops.tsv')
        args.phi_base       = d(args.phi_base,       'phi_base/phi_base.tsv')
        args.dbcan_tc       = d(args.dbcan_tc,       'dbcan/diamond.out.tc')
    return args


def main():
    args = parse_args()
    annotations = defaultdict(new_entry)
    mrna_to_gene = {}

    # ── ID map (old pre-rename ID → new ID) ─────────────────────────────────
    # Tool outputs (eggNOG, IPS6, dbCAN, ...) were run against the protein
    # FASTA under the OLD, pre-rename_gff_ids.py ID namespace. Translate
    # their query IDs into the current GFF's namespace before resolving.
    id_translate = {}
    if args.id_map and not is_absent(args.id_map):
        with open(args.id_map) as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 2 and c[0] != c[1]:
                    id_translate[c[0]] = c[1]

    # ── 0. Parse GFF3 topology ──────────────────────────────────────────────
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            ftype, attrs = cols[2], cols[8]
            id_m = re.search(r'ID=([^;]+)', attrs)
            par_m = re.search(r'Parent=([^;]+)', attrs)
            if ftype == 'gene' and id_m:
                _ = annotations[id_m.group(1)]
            elif ftype in ('mRNA', 'ncRNA', 'transcript') and id_m and par_m:
                mid, gid = id_m.group(1), par_m.group(1)
                mrna_to_gene[mid] = gid
                annotations[gid]['mrna_ids'].add(mid)

    def resolve(qid):
        qid = id_translate.get(qid, qid)
        if qid in mrna_to_gene:
            return mrna_to_gene[qid]
        if qid in annotations:
            return qid
        stripped = re.sub(r'[-_.][Tt]\d+$|\.mrna\d+$|\.t\d+$', '', qid)
        stripped = id_translate.get(stripped, stripped)
        if stripped in mrna_to_gene:
            return mrna_to_gene[stripped]
        return stripped if stripped in annotations else qid

    # ── 1. eggNOG-mapper annotations ────────────────────────────────────────
    # emapper v3/e7 22-col header (eggnog-mapper 3.0.0-beta6+): #query
    # seed_ortholog evalue score eggNOG_OGs tax_ceiling farthest_donor_lineage
    # COG_category Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module
    # KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs
    # annotation_confidence.
    # v3 dropped the old v2/beta5 `Description` column outright (upstream:
    # "Preferred_name is retained and is the useful human-readable gene label
    # for e7") — there is no replacement source for it, so nog_desc/
    # eggNOG_description is left unpopulated from eggNOG going forward.
    if not is_absent(args.eggnog):
        with open(args.eggnog) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 22:
                    continue
                a = annotations[resolve(c[0])]
                # Neither eggNOG field drives product/gene_symbol directly —
                # priority is KofamScan first, then eggNOG's KO (via KEGG
                # REST, in step 3b), never eggNOG's own Preferred_name. It's
                # gated behind eggNOG's own internal ortholog-confidence
                # threshold in an inconsistent way that produced real
                # mislabeling (e.g. a gene eggNOG's own Description called
                # "Enoyl-(Acyl carrier protein) reductase" getting symbol
                # "SOU1" from a broad, promiscuous KO match) — excluded
                # entirely.
                if c[9]  and c[9]  != '-': a['go'].update(g.strip() for g in c[9].split(',') if g.strip())
                if c[10] and c[10] != '-': a['ec'].update(e.strip() for e in c[10].split(',') if e.strip())
                if c[11] and c[11] != '-':
                    kos = [k.strip().removeprefix('ko:') for k in c[11].split(',') if k.strip()]
                    a['kegg_ko'].update(kos)
                    if kos and not a['ko_primary']:
                        a['ko_primary'] = kos[0]   # lowest priority — KofamScan overwrites this below
                if c[12] and c[12] != '-': a['kegg_path'].update(p.strip() for p in c[12].split(',') if p.strip())
                if c[4]  and c[4]  != '-': a['nog_og'] = a['nog_og'] or c[4]
                if c[7]  and c[7]  != '-': a['cog'] = a['cog'] or c[7]

    # ── 2. InterProScan6 TSV ────────────────────────────────────────────────
    # Cols: 0=prot 3=analysis 4=sig_acc 5=sig_desc 6=start 7=stop 8=score
    #       11=ipr_acc 13=GO 14=pathways
    # TMbed descriptor lives in col4 (sig_acc), not col5 — col5 is always '-'
    # for TMbed rows. Values: TMhelix_in-to-out / TMhelix_out-to-in /
    # TMbeta_in-to-out / TMbeta_out-to-in / Signal_peptide.
    if not is_absent(args.interpro):
        tm_tmbed = defaultdict(int)
        tm_phobius = defaultdict(int)
        sp_phobius = set()
        sp_tmbed = set()
        with open(args.interpro) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 5:
                    continue
                prot_id = c[0]
                a = annotations[resolve(prot_id)]
                anal, sacc = c[3], c[4]
                if len(c) > 11 and c[11] and c[11] != '-':
                    a['interpro'].add(c[11])
                if len(c) > 13 and c[13] and c[13] != '-':
                    # InterPro's GO column carries a source suffix per term,
                    # e.g. "GO:0005741(InterPro)|GO:0005741(PANTHER)" — the
                    # same term from two source tools. Strip it so identical
                    # GO IDs collapse into one set entry instead of duplicating.
                    a['go'].update(g.strip().split('(')[0] for g in c[13].split('|') if g.strip().startswith('GO:'))
                if anal == 'Pfam' and sacc and sacc.startswith('PF'):
                    a['pfam'].add(sacc)
                if anal == 'PANTHER' and sacc:
                    a['panther'] = a['panther'] or sacc
                if anal.lower().startswith('signalp'):
                    a['signalp'] = True
                if anal == 'TMbed':
                    if sacc.startswith('TMhelix') or sacc.startswith('TMbeta'):
                        tm_tmbed[prot_id] += 1
                    elif sacc == 'Signal_peptide':
                        sp_tmbed.add(prot_id)
                if anal == 'Phobius':
                    if sacc == 'TRANSMEMBRANE':
                        tm_phobius[prot_id] += 1
                    elif 'SIGNAL' in sacc:
                        sp_phobius.add(prot_id)
        for prot_id, cnt in tm_tmbed.items():
            gid = resolve(prot_id)
            annotations[gid]['tm_tmbed'] = max(annotations[gid]['tm_tmbed'], cnt)
        for prot_id, cnt in tm_phobius.items():
            gid = resolve(prot_id)
            annotations[gid]['tm_phobius'] = max(annotations[gid]['tm_phobius'], cnt)
        for prot_id in sp_phobius:
            annotations[resolve(prot_id)]['sp_phobius'] = True
        for prot_id in sp_tmbed:
            annotations[resolve(prot_id)]['signalp'] = True

    # ── 3. KofamScan ────────────────────────────────────────────────────────
    # KofamScan has priority over eggNOG for product/ko_primary: its own KO
    # overwrites any eggNOG-set placeholder, but only the FIRST KofamScan hit
    # per gene wins if a gene has more than one above-threshold match.
    if not is_absent(args.kofamscan):
        kofam_primary_set = set()
        with open(args.kofamscan) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 7 or c[0] != '*':
                    continue
                a = annotations[resolve(c[1])]
                gid_key = resolve(c[1])
                if c[2]:
                    a['kegg_ko'].add(c[2])
                    if gid_key not in kofam_primary_set:
                        a['ko_primary'] = c[2]
                # KofamScan's detail-tsv wraps the KO definition in double
                # quotes, e.g. "alpha-L-fucosidase 2 [EC:3.2.1.51]". Strip
                # them FIRST — otherwise the trailing " sits after the ], so
                # ec_tail.rstrip(']') can't remove the bracket and the EC
                # comes out malformed as e.g. 3.2.1.51]".
                desc = c[6].strip().strip('"').strip()
                if '[EC:' in desc:
                    prod, _, ec_tail = desc.partition('[EC:')
                    for ec in ec_tail.rstrip(']').split():
                        a['ec'].add(ec)
                    desc = prod.strip()
                if desc and gid_key not in kofam_primary_set:
                    a['product'] = desc
                kofam_primary_set.add(gid_key)

    # ── 3b. KEGG REST enrichment (bin/fetch_kegg_ko_info.py) ────────────────
    # KO-keyed lookup (gene_symbols, ko_name), covering KOs from BOTH
    # KofamScan and eggNOG. Driven by a['ko_primary'] — KofamScan's KO if it
    # hit this gene, else eggNOG's KO — never by eggNOG's own Preferred_name/
    # Description text directly (see step 1's comment for why: eggNOG's own
    # confidence gating on Preferred_name proved inconsistent and produced
    # real mislabeling). gene_symbol always comes from here when a KO exists;
    # product only when KofamScan didn't already supply one.
    if not is_absent(args.ko_info):
        ko_lookup = {}
        with open(args.ko_info) as fh:
            next(fh, None)
            for line in fh:
                if not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 4:
                    continue
                ko_lookup[c[0]] = {'symbols': c[1], 'name': c[2]}
        for a in annotations.values():
            if not a['ko_primary']:
                continue
            info = ko_lookup.get(a['ko_primary'])
            if not info:
                continue
            if info['symbols'] and not a['gene_symbol']:
                a['gene_symbol'] = info['symbols'].split(',')[0].strip()
            if info['name'] and not a['product']:
                # KEGG's NAME field embeds EC numbers inline, e.g.
                # "amidase [EC:3.5.1.4]" — strip it out (same as the
                # KofamScan path) since it's already captured separately
                # in a['ec'], to avoid duplicating it in the product text.
                name = info['name'].strip()
                if '[EC:' in name:
                    prod, _, ec_tail = name.partition('[EC:')
                    for ec in ec_tail.rstrip(']').split():
                        a['ec'].add(ec)
                    name = prod.strip()
                a['product'] = name

    # ── 3c. gene2product curated name lookup ────────────────────────────────
    # Runs after gene_symbol is actually populated (step 3b) — a curated
    # override keyed on gene_symbol can't match anything any earlier.
    if args.gene2product and not is_absent(args.gene2product):
        _g2p = {}
        with open(args.gene2product) as fh:
            for ln in fh:
                if ln.startswith('#') or not ln.strip():
                    continue
                pcs = ln.rstrip('\n').split('\t')
                if len(pcs) >= 2:
                    _g2p[pcs[0].upper()] = pcs[1]
        if _g2p:
            for a in annotations.values():
                sym = a['gene_symbol']
                if not sym:
                    continue
                hit = _g2p.get(sym.upper())
                if hit and (not a['product'] or a['product'] == sym):
                    a['product'] = hit

    # ── 4. Infernal tblout ──────────────────────────────────────────────────
    # --fmt 2 --tblout cols: 0=idx 1=target_name 2=accession(RFxxxxx)
    # 3=query_name(our gene/mRNA ID) 4=query_accession ...
    if not is_absent(args.infernal):
        with open(args.infernal) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.split()
                if len(c) >= 4:
                    annotations[resolve(c[3])]['rfam'].add(c[2])

    # ── 5. EffectorP3 ───────────────────────────────────────────────────────
    # Cols: 0=Identifier 1=Cytoplasmic_effector 2=Apoplastic_effector
    #       3=Non-effector 4=Prediction. "Y"-cells carry a probability suffix
    #       e.g. "Y (0.92)" — must check startswith, not exact-match.
    if not is_absent(args.effectorp3):
        with open(args.effectorp3) as fh:
            for line in fh:
                line = line.rstrip()
                if not line or line.startswith('#') or line.startswith('Identifier'):
                    continue
                c = line.split('\t')
                if len(c) < 3:
                    continue
                if c[1].strip().startswith('Y') or c[2].strip().startswith('Y'):
                    a = annotations[resolve(c[0].split()[0])]
                    a['effector'] = c[4].strip() if len(c) > 4 else 'Effector'

    # ── 6. dbCAN overview ───────────────────────────────────────────────────
    if not is_absent(args.dbcan_overview):
        with open(args.dbcan_overview) as fh:
            next(fh, None)
            for line in fh:
                if not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 7:
                    continue
                a = annotations[resolve(c[0])]
                n_tools = int(c[5]) if c[5].isdigit() else 0
                family = c[6].strip()
                if not family or family in ('-', 'N', '0'):
                    for fc in c[2:5]:
                        fc = fc.strip()
                        if fc and fc not in ('-', 'N', '0'):
                            family = fc.split('(')[0].strip()
                            break
                if family and family not in ('-', 'N', '0'):
                    a['cazyme'] = family
                    a['caz_tools'] = n_tools

    # ── 7. dbCAN CGC / PUL clusters ─────────────────────────────────────────
    if not is_absent(args.dbcan_cgc):
        with open(args.dbcan_cgc) as fh:
            next(fh, None)
            for line in fh:
                if not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) >= 4 and c[3]:
                    annotations[resolve(c[3])]['pul'] = c[0]

    # ── 8. MEROPS blastp ────────────────────────────────────────────────────
    if not is_absent(args.merops):
        with open(args.merops) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 7:
                    continue
                a = annotations[resolve(c[0])]
                if a['mer_hit']:
                    continue
                a['mer_hit'] = c[1]
                fam_m = re.search(r'\[([A-Z]\d+\.\d+)\]', c[6])
                if fam_m:
                    a['mer_fam'] = fam_m.group(1)

    # ── 8b. dbCAN diamond.out.tc — TCDB transporter classification ──────────
    # Cols: 0=TCDB_ID 1=TCDB_len 2=query(gene) 3=query_len 4=evalue
    #       5=TCDB_start 6=TCDB_end 7=qstart 8=qend 9=coverage 10=Database
    if not is_absent(args.dbcan_tc):
        with open(args.dbcan_tc) as fh:
            for line in fh:
                if line.startswith('TCDB ID') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 3:
                    continue
                a = annotations[resolve(c[2])]
                if not a['tc']:
                    a['tc'] = c[0]

    # ── 9. PHI-base blastp ───────────────────────────────────────────────────
    if not is_absent(args.phi_base):
        with open(args.phi_base) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 2:
                    continue
                a = annotations[resolve(c[0])]
                if a['phi_acc']:
                    continue
                header = c[6] if len(c) > 6 else c[1]
                phi_p = header.split('#')
                # PHI-base joins multiple values for one hit with '__'
                # (e.g. "PHI:1816__PHI:12191" for a protein with two PHI
                # entries) — normalize to a comma-separated list.
                if len(phi_p) >= 2 and phi_p[1].startswith('PHI:'):
                    a['phi_acc'] = phi_p[1].replace('__', ',')
                if len(phi_p) >= 3:
                    a['phi_gene'] = phi_p[2].replace('__', ',')
                if len(phi_p) >= 6:
                    a['phi_pheno'] = phi_p[5].replace('__', ',')

    # ── 10. AntiSMASH BGC annotations (from combined GBK) ───────────────────
    if args.antismash_gbk and not is_absent(args.antismash_gbk):
        b = {'l': None, 'k': '', 'f': [], 'd': set(), 'g': [], 'q': None, 'v': '', 's': False}

        def commit_q(ft):
            if ft != 'CDS' or not b['q']:
                b['q'] = None; b['v'] = ''; b['s'] = False; return
            v, q = b['v'], b['q']
            if q == 'locus_tag':
                # antiSMASH appends an 8-char hex hash to locus_tag for CDS
                # features duplicated across overlapping cluster regions in
                # its combined GBK (e.g. "SS02_016147_82dc74fc") — strip it
                # so it resolves back to the real gene ID.
                b['l'] = re.sub(r'_[0-9a-f]{8}$', '', v)
            elif q == 'gene_kind':
                b['k'] = v
            elif q == 'gene_functions':
                b['f'].append(v)
                sm = re.search(r'\(smcogs\)\s+(SMCOG\d+:[^(]+)', v)
                if sm: b['g'].append(sm.group(1).strip())
            elif q == 'note':
                sm = re.search(r'smCOG[:\s]\s*(SMCOG\d+:[^(]+)', v, re.IGNORECASE)
                if sm: b['g'].append(sm.group(1).strip())
            elif q == 'sec_met_domain':
                dom = v.split(' (')[0].strip()
                if dom: b['d'].add(dom)
            b['q'] = None; b['v'] = ''; b['s'] = False

        def store():
            if b['l'] and (b['k'] or b['f'] or b['g']):
                bt = ''
                for fn in b['f']:
                    bm = re.search(r'rule-based-clusters\)\s+(\S+?):', fn)
                    if bm: bt = (bt + '/' if bt else '') + bm.group(1)
                a = annotations[resolve(b['l'])]
                if bt: a['bgc_type'] = bt
                if b['k']: a['bgc_role'] = b['k']
                a['bgc_domains'].update(b['d'])
                if b['g']: a['smcog'] = b['g'][0]
            b.update({'l': None, 'k': '', 'f': [], 'd': set(), 'g': [], 'q': None, 'v': '', 's': False})

        bft = None
        with open(args.antismash_gbk) as afh:
            for ar in afh:
                al = ar.rstrip('\n')
                if al.startswith('ORIGIN') or al == '//':
                    if bft == 'CDS': commit_q(bft); store()
                    bft = None
                elif re.match(r'^     \S', al):
                    if bft == 'CDS': commit_q(bft); store()
                    bft = al.split()[0]
                elif al.startswith('                     /') and not b['s']:
                    commit_q(bft)
                    rest = al[22:]
                    if '=' in rest:
                        b['q'], _, vv = rest.partition('=')
                        vv = vv.lstrip('"')
                        if vv.endswith('"'): b['v'] = vv[:-1]
                        else: b['v'] = vv; b['s'] = True
                    else:
                        b['q'] = rest.strip(); b['v'] = ''
                elif (b['s'] and al.startswith('                     ')
                      and not al.startswith('                     /')):
                    ct = al.strip()
                    if ct.endswith('"'): b['v'] += ' ' + ct[:-1]; b['s'] = False
                    else: b['v'] += ' ' + ct

    # ── GO term propagation via go-basic.obo ────────────────────────────────
    if args.go_obo and not is_absent(args.go_obo):
        par, alt = {}, {}
        cur_id, cur_par, cur_alt, cur_obs = None, set(), [], False
        with open(args.go_obo) as fh:
            for raw in fh:
                ln = raw.rstrip()
                if ln == '[Term]':
                    if cur_id and not cur_obs:
                        par[cur_id] = cur_par
                        for a_ in cur_alt: alt[a_] = cur_id
                    cur_id, cur_par, cur_alt, cur_obs = None, set(), [], False
                elif ln.startswith('id: GO:'):
                    cur_id = ln[4:]
                elif ln.startswith('alt_id: GO:'):
                    cur_alt.append(ln[8:])
                elif ln == 'is_obsolete: true':
                    cur_obs = True
                elif not cur_obs and cur_id:
                    if ln.startswith('is_a: '):
                        cur_par.add(ln[6:].split()[0])
                    elif ln.startswith('relationship: part_of '):
                        pp = ln.split()[2]
                        if pp.startswith('GO:'): cur_par.add(pp)
        if cur_id and not cur_obs:
            par[cur_id] = cur_par
            for a_ in cur_alt: alt[a_] = cur_id

        anc = {}
        def ancestors(goid):
            if goid in anc: return anc[goid]
            result, stack = set(), list(par.get(goid, set()))
            while stack:
                p = stack.pop()
                if p not in result:
                    result.add(p)
                    stack.extend(par.get(p, set()) - result)
            anc[goid] = result
            return result

        for a in annotations.values():
            if a['go']:
                prop = set(a['go'])
                for t in list(a['go']):
                    can = alt.get(t, t)
                    prop.update(ancestors(can))
                    if can != t: prop.add(can)
                a['go_propagated'] = prop
            else:
                a['go_propagated'] = set()
    else:
        for a in annotations.values():
            a['go_propagated'] = set(a['go'])

    # ── Write functional_annotation.tsv ─────────────────────────────────────
    COLS = [
        'gene_id', 'mrna_ids', 'product', 'gene_symbol',
        'GO_terms', 'GO_propagated_count', 'EC_numbers', 'KEGG_KO', 'KEGG_pathways',
        'COG_category', 'eggNOG_OG', 'eggNOG_description',
        'InterPro_accessions', 'Pfam_domains', 'PANTHER',
        'Secreted', 'SignalP', 'TM_helices_TMbed', 'TM_helices_Phobius', 'SP_Phobius',
        'EffectorP_class',
        'CAZyme_family', 'CAZyme_tools_agree', 'PUL_cluster',
        'MEROPS_hit', 'MEROPS_family', 'TC_id',
        'PHIbase_accession', 'PHIbase_gene', 'PHIbase_phenotype',
        'Rfam_accessions',
        'BGC_cluster_type', 'BGC_gene_role', 'BGC_domains', 'SMCOG',
    ]
    with open(args.out_tsv, 'w') as out:
        out.write('\t'.join(COLS) + '\n')
        for gid in sorted(annotations):
            a = annotations[gid]
            row = [
                gid, ','.join(sorted(a['mrna_ids'])), ncbi_product(a['product'], a['gene_symbol']), a['gene_symbol'],
                '|'.join(sorted(a['go'])), str(len(a['go_propagated'])),
                ','.join(sorted(a['ec'])), ','.join(sorted(a['kegg_ko'])), '|'.join(sorted(a['kegg_path'])),
                a['cog'], a['nog_og'], a['nog_desc'],
                '|'.join(sorted(a['interpro'])), '|'.join(sorted(a['pfam'])), a['panther'],
                'Y' if is_secreted(a) else '', 'Y' if a['signalp'] else '',
                str(a['tm_tmbed']), str(a['tm_phobius']), 'Y' if a['sp_phobius'] else '',
                a['effector'], a['cazyme'], str(a['caz_tools']), a['pul'],
                a['mer_hit'], a['mer_fam'], a['tc'],
                a['phi_acc'], a['phi_gene'], a['phi_pheno'],
                '|'.join(sorted(a['rfam'])),
                a['bgc_type'], a['bgc_role'], '|'.join(sorted(a['bgc_domains'])), a['smcog'],
            ]
            out.write('\t'.join(row) + '\n')

    # ── Write annotated GFF3 ─────────────────────────────────────────────────
    with open(args.gff) as fin, open(args.out_gff, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                fout.write(line)
                continue
            ftype, attrs = cols[2], cols[8]
            id_m = re.search(r'ID=([^;]+)', attrs)
            par_m = re.search(r'Parent=([^;]+)', attrs)

            is_transcript = ftype in ('mRNA', 'ncRNA', 'transcript')

            if ftype == 'gene':
                if not id_m:
                    fout.write(line); continue
                gid = id_m.group(1)
            elif is_transcript:
                if not id_m or not par_m:
                    fout.write(line); continue
                gid = mrna_to_gene.get(id_m.group(1), par_m.group(1))
            elif ftype == 'CDS':
                if not par_m:
                    fout.write(line); continue
                parent_id = par_m.group(1).split(',')[0]
                gid = mrna_to_gene.get(parent_id, parent_id)
            else:
                fout.write(line); continue

            if gid not in annotations:
                fout.write(line); continue

            a = annotations[gid]
            if ftype == 'ncRNA' and not a['product']:
                product = 'hypothetical lncRNA'
            else:
                product = ncbi_product(a['product'], a['gene_symbol'])
            symbol = a['gene_symbol']
            go_set = a['go_propagated'] or a['go']
            secreted = is_secreted(a)
            extras = []

            # NCBI-approved /db_xref databases only (InterPro, PFAM, RFAM).
            # KEGG and COG are not on NCBI's approved db_xref list — they go
            # into Note instead. GO already gets its own Ontology_term qualifier.
            dbx = (['InterPro:' + i for i in sorted(a['interpro'])]
                   + ['PFAM:' + p for p in sorted(a['pfam'])]
                   + ['RFAM:' + r for r in sorted(a['rfam'])])

            if ftype == 'gene':
                if symbol:
                    extras.append(('Name', gff_val(symbol)))
            elif is_transcript:
                # product, EC_number, and Dbxref are transcript-level only
                # (mRNA/ncRNA) — not on CDS. No gene= qualifier — the gene
                # symbol lives on the gene feature's Name= only.
                extras.append(('product', gff_val(product)))
                if go_set:
                    extras.append(('Ontology_term', ','.join(sorted(go_set))))
                if a['ec']:
                    extras.append(('EC_number', ','.join(sorted(a['ec']))))
                if dbx:
                    extras.append(('Dbxref', ','.join(dbx)))

            # Note is transcript-level only (mRNA/ncRNA) — not on gene or CDS.
            # Entries are comma-separated, unescaped (no percent-encoding) —
            # none of the constructed fragments below contain a literal ';'
            # or '=', so no GFF3 attribute-grammar character needs escaping.
            # Multi-value sub-fields (KEGG, Rfam) use '|' internally so they
            # can't be mistaken for separate Note entries at the comma level.
            if is_transcript:
                notes = []
                if secreted: notes.append('Secreted')
                if a['tm_tmbed'] > 0: notes.append(f"Tmbed:TM{a['tm_tmbed']}")
                if a['tm_phobius'] > 0: notes.append(f"Phobius:TM{a['tm_phobius']}")
                if a['effector']: notes.append(a['effector'])
                if a['cazyme']: notes.append(f"CAZyme:{a['cazyme']}")
                if a['pul']: notes.append(f"Polysaccharide utilization locus {a['pul']}")
                if a['mer_hit']: notes.append(f"MEROPS:{a['mer_hit']}")
                if a['tc']: notes.append(f"TC:{a['tc']}")
                if a['phi_acc']: notes.append(a['phi_acc'])
                if a['cog']: notes.append(f"COG:{a['cog']}")
                notes.extend(f"KEGG:{ko}" for ko in sorted(a['kegg_ko']))
                if a['bgc_type']:
                    notes.append(f"BGC:{a['bgc_type']}")
                elif a['bgc_role'] and a['bgc_role'] != 'other':
                    notes.append(f"BGC:{a['bgc_role']}")
                if a['smcog']:
                    smcog_m = re.search(r'SMCOG(\d+)', a['smcog'])
                    notes.append(f"SMCOG:{smcog_m.group(1)}" if smcog_m else f"SMCOG:{a['smcog']}")
                if notes:
                    extras.append(('Note', ','.join(notes)))

            cols[8] = merge_attrs(attrs, extras)
            fout.write('\t'.join(cols) + '\n')

    # ── Write stats ───────────────────────────────────────────────────────────
    total = len(annotations)
    def pct(n):
        return '%.1f%%' % (n * 100.0 / total) if total else '0.0%'

    rows = [
        ('Total genes', total),
        ('With product name', sum(1 for a in annotations.values() if a['product'])),
        ('With GO terms (direct)', sum(1 for a in annotations.values() if a['go'])),
        ('With GO terms (propagated)', sum(1 for a in annotations.values() if a['go_propagated'])),
        ('With KEGG KO', sum(1 for a in annotations.values() if a['kegg_ko'])),
        ('With EC numbers', sum(1 for a in annotations.values() if a['ec'])),
        ('With InterPro domains', sum(1 for a in annotations.values() if a['interpro'])),
        ('With Pfam domains', sum(1 for a in annotations.values() if a['pfam'])),
        ('Secreted proteins', sum(1 for a in annotations.values() if is_secreted(a))),
        ('Signal peptide (SignalP)', sum(1 for a in annotations.values() if a['signalp'])),
        ('TM proteins (TMbed)', sum(1 for a in annotations.values() if a['tm_tmbed'] > 0)),
        ('TM proteins (Phobius)', sum(1 for a in annotations.values() if a['tm_phobius'] > 0)),
        ('Predicted effectors', sum(1 for a in annotations.values() if a['effector'])),
        ('CAZymes (dbCAN)', sum(1 for a in annotations.values() if a['cazyme'])),
        ('Genes in PUL clusters', sum(1 for a in annotations.values() if a['pul'])),
        ('MEROPS peptidases', sum(1 for a in annotations.values() if a['mer_hit'])),
        ('Transporters (TCDB)', sum(1 for a in annotations.values() if a['tc'])),
        ('PHI-base hits', sum(1 for a in annotations.values() if a['phi_acc'])),
        ('ncRNA genes (Rfam)', sum(1 for a in annotations.values() if a['rfam'])),
        ('Genes in BGC clusters', sum(1 for a in annotations.values() if a['bgc_type'] or a['bgc_role'])),
        ('Genes with smCOG hit', sum(1 for a in annotations.values() if a['smcog'])),
    ]
    with open(args.out_stats, 'w') as out:
        out.write('=== bagRNA Functional Annotation Summary ===\n\n')
        for label, n in rows:
            out.write('%-32s  %6d  (%s)\n' % (label, n, pct(n)))

    print(f"Wrote {args.out_gff}, {args.out_tsv}, {args.out_stats}")
    print(f"Total genes annotated: {total}")


if __name__ == '__main__':
    main()
