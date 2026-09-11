#!/usr/bin/env python3
"""
Build a summary table (.tbl), a self-contained HTML report, and a PDF
rendering of that same HTML (via headless Chromium) over a completed bagRNA
annotation run.

Every input past --genome-fasta/--final-gff/--proteins-faa is optional and
follows this pipeline's NO_FILE sentinel convention (an absent/empty path
means that pipeline stage was skipped) — each corresponding report section
is simply left out rather than failing the whole report.
"""
import argparse
import html
import json
import os
import shutil
import subprocess

COG_NAMES = {
    "J": "Translation, ribosomal structure", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics", "D": "Cell cycle control, cell division",
    "Y": "Nuclear structure", "V": "Defense mechanisms", "T": "Signal transduction",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "Z": "Cytoskeleton", "W": "Extracellular structures",
    "U": "Intracellular trafficking and secretion", "O": "PTM, protein turnover, chaperones",
    "C": "Energy production and conversion", "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism", "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism", "Q": "Secondary metabolite biosynthesis",
    "R": "General function prediction only", "S": "Function unknown",
    "X": "Mobilome: prophages, transposons",
}


def is_absent(path):
    return not path or not os.path.exists(path) or os.path.getsize(path) == 0


def trunc(s, n=70):
    return s if len(s) <= n else s[: n - 1] + "…"


def esc(s):
    return html.escape(str(s), quote=True)


def fmt(n):
    return f"{n:,}" if isinstance(n, (int, float)) else str(n)


def compact_bp(n):
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f} Mb"
    if n >= 1_000:
        return f"{n / 1_000:.1f} kb"
    return f"{n} bp"


# --------------------------------------------------------------- parsers ---
def genome_stats(path):
    seqs, cur = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append(len("".join(cur)))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append(len("".join(cur)))
    seqs.sort(reverse=True)
    total = sum(seqs)
    csum, n50 = 0, None
    for L in seqs:
        csum += L
        if csum >= total / 2 and n50 is None:
            n50 = L
    return {"contigs": len(seqs), "total_bp": total, "n50": n50, "largest": seqs[0] if seqs else 0}


def busco_stats(path):
    stats = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("## lineage:"):
                stats["lineage"] = line.split(":", 1)[1].strip()
            elif ":" in line and "," in line:
                key, rest = line.split(":", 1)
                pct, _, n = rest.partition(",")
                stats[key.strip()] = float(pct.strip().rstrip("%"))
            elif line.startswith("N:"):
                stats["N"] = int(line.split(":", 1)[1].strip())
    return stats if {"S", "D", "F", "M"} <= stats.keys() else None


def feature_counts(gff_path):
    counts = {}
    gene_lengths, ncrna_lengths = [], []
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            ftype = cols[2]
            counts[ftype] = counts.get(ftype, 0) + 1
            try:
                length = int(cols[4]) - int(cols[3]) + 1
            except ValueError:
                continue
            if ftype == "gene":
                gene_lengths.append(length)
            elif ftype == "ncRNA":
                ncrna_lengths.append(length)
    return counts, gene_lengths, ncrna_lengths


def protein_count(faa_path):
    with open(faa_path) as f:
        return sum(1 for l in f if l.startswith(">"))


def functional_stats(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("==="):
                continue
            parts = line.rsplit(None, 2)
            if len(parts) == 3:
                label, count, pct = parts
                try:
                    rows.append((label.strip(), int(count), float(pct.strip("()%"))))
                except ValueError:
                    continue
    return rows


def cazyme_classes(overview_path):
    counts = {}
    n_genes = 0
    with open(overview_path) as f:
        next(f, None)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            n_genes += 1
            val = cols[2] if cols[2] != "-" else cols[3]
            for fam in val.split("+"):
                fam = fam.split("(")[0].split("_")[0]
                cls = "".join(ch for ch in fam if ch.isalpha())
                if cls:
                    counts[cls] = counts.get(cls, 0) + 1
    return counts, n_genes


def antismash_bgc_types(json_path):
    counts = {}
    n_regions = 0
    with open(json_path) as f:
        d = json.load(f)
    for r in d.get("records", []):
        for feat in r.get("features", []):
            if feat.get("type") == "region":
                n_regions += 1
                for p in feat.get("qualifiers", {}).get("product", []):
                    counts[p] = counts.get(p, 0) + 1
    return counts, n_regions


def phi_base_phenotypes(tsv_path, top_n=8):
    counts = {}
    with open(tsv_path) as f:
        for line in f:
            pheno = line.rstrip("\n").split("#")[-1]
            if pheno:
                counts[pheno] = counts.get(pheno, 0) + 1
    top = sorted(counts.items(), key=lambda x: -x[1])[:top_n]
    return top, sum(counts.values())


def cog_categories(eggnog_path, cog_def_path):
    cog_id_to_cat = {}
    if not is_absent(cog_def_path):
        with open(cog_def_path, encoding="latin-1") as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    cog_id_to_cat[cols[0]] = cols[1]

    counts = {}
    n_unmapped = n_none = 0
    with open(eggnog_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            val = cols[7]
            if val in ("-", ""):
                n_none += 1
                continue
            if val.startswith("COG"):
                cats = cog_id_to_cat.get(val)
                if not cats:
                    n_unmapped += 1
                    continue
            else:
                cats = val
            for c in cats:
                if c in COG_NAMES:
                    counts[c] = counts.get(c, 0) + 1
    return counts, n_unmapped, n_none


def kegg_modules(json_path, min_total=4, top_n=15):
    with open(json_path) as f:
        modules = json.load(f)
    sized = [m for m in modules if m["total"] >= min_total]
    sized.sort(key=lambda x: (-x["pct"], -x["total"]))
    return sized, sized[:top_n]


# --------------------------------------------------------- HTML building ---
CSS = """
.br-root {
  color-scheme: light;
  --surface-1:      #fcfcfb;
  --page-plane:     #f9f9f7;
  --text-primary:   #0b0b0b;
  --text-secondary: #52514e;
  --text-muted:     #898781;
  --grid:           #e1e0d9;
  --border:         rgba(11,11,11,0.10);
  --series-1:       #2a78d6;
  --series-1-track: #dbe8f8;
  /* BUSCO composition (single/duplicated/fragmented/missing): categorical
     slots 1-4 in their fixed documented order — validated as a set
     (adjacent-pair CVD/contrast) via the dataviz skill's validator, rather
     than status colors repurposed as four arbitrary series (which fails:
     good-green and duplicated-aqua sit under Delta E 15 apart). */
  --busco-1:        #2a78d6;
  --busco-2:        #eb6834;
  --busco-3:        #1baf7a;
  --busco-4:        #eda100;
  --status-good:    #0ca30c;
}
@media (prefers-color-scheme: dark) {
  .br-root {
    color-scheme: dark;
    --surface-1:      #1a1a19;
    --page-plane:     #0d0d0d;
    --text-primary:   #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted:     #898781;
    --grid:           #2c2c2a;
    --border:         rgba(255,255,255,0.10);
    --series-1:       #3987e5;
    --series-1-track: #23324a;
    --busco-1:        #3987e5;
    --busco-2:        #d95926;
    --busco-3:        #199e70;
    --busco-4:        #c98500;
    --status-good:    #0ca30c;
  }
}
.br-root * { box-sizing: border-box; }
.br-root {
  background: var(--page-plane);
  color: var(--text-primary);
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  padding: 32px 20px 56px;
}
.br-wrap { max-width: 940px; margin: 0 auto; }
.br-header { margin-bottom: 28px; }
.br-header h1 { font-size: 26px; font-weight: 700; margin: 0 0 4px; }
.br-header .br-species { font-size: 16px; font-style: italic; color: var(--text-secondary); margin: 0 0 2px; }
.br-header .br-sub { font-size: 12px; color: var(--text-muted); margin: 0; }

.br-kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(168px, 1fr)); gap: 12px; margin-bottom: 28px; }
.br-tile { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 14px 16px; }
.br-tile .br-tile-label { font-size: 11.5px; color: var(--text-secondary); margin-bottom: 6px; }
.br-tile .br-tile-value { font-size: 20px; font-weight: 600; color: var(--text-primary); white-space: nowrap; }

.br-card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 20px 22px; margin-bottom: 20px; }
.br-card h2 { font-size: 15px; font-weight: 600; margin: 0 0 2px; }
.br-card .br-card-sub { font-size: 11.5px; color: var(--text-muted); margin: 0 0 16px; }
.br-card .br-footnote { font-size: 10.5px; color: var(--text-muted); margin-top: 12px; line-height: 1.5; }

.br-bar-row { display: grid; grid-template-columns: minmax(120px, 230px) 1fr auto; gap: 12px; align-items: center; padding: 3px 0; }
.br-bar-label { font-size: 12px; color: var(--text-secondary); overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.br-bar-track { background: var(--grid); border-radius: 4px; height: 14px; overflow: hidden; }
.br-bar-fill { height: 100%; border-radius: 0 4px 4px 0; background: var(--series-1); }
.br-bar-value { font-size: 12px; color: var(--text-primary); font-variant-numeric: tabular-nums; white-space: nowrap; min-width: 70px; text-align: right; }

.br-stack { display: flex; height: 22px; border-radius: 5px; overflow: hidden; gap: 2px; background: var(--page-plane); }
.br-stack-seg { height: 100%; }
.br-legend { display: flex; flex-wrap: wrap; gap: 14px; margin-top: 12px; }
.br-legend-item { display: flex; align-items: center; gap: 6px; font-size: 12px; color: var(--text-secondary); }
.br-legend-dot { width: 10px; height: 10px; border-radius: 3px; flex: none; }

.br-footer { text-align: center; font-size: 11px; color: var(--text-muted); margin-top: 8px; }
"""


def bar_section(title, subtitle, items, unit="", footnote=None):
    """items: list of (label, value, display_text)."""
    if not items:
        return ""
    max_v = max(v for _, v, _ in items) or 1
    rows = []
    for label, value, disp in items:
        pct = max(2, round(100 * value / max_v))
        rows.append(f"""
        <div class="br-bar-row">
          <div class="br-bar-label" title="{esc(label)}">{esc(label)}</div>
          <div class="br-bar-track"><div class="br-bar-fill" style="width:{pct}%"></div></div>
          <div class="br-bar-value">{esc(disp)}</div>
        </div>""")
    foot = f'<p class="br-footnote">{esc(footnote)}</p>' if footnote else ""
    return f"""
    <div class="br-card">
      <h2>{esc(title)}</h2>
      <p class="br-card-sub">{esc(subtitle)}</p>
      {"".join(rows)}
      {foot}
    </div>"""


def kegg_module_section(title, subtitle, modules, footnote):
    if not modules:
        return ""
    rows = []
    for m in modules:
        pct = m["pct"]
        is_complete = pct >= 100
        color = "var(--status-good)" if is_complete else "var(--series-1)"
        label = f"{m['module_id']} {trunc(m['module_name'])}"
        disp = f"{pct}% ({m['found']}/{m['total']} steps)" + (" ✓" if is_complete else "")
        rows.append(f"""
        <div class="br-bar-row">
          <div class="br-bar-label" title="{esc(label)}">{esc(label)}</div>
          <div class="br-bar-track"><div class="br-bar-fill" style="width:{max(2, round(pct))}%; background:{color}"></div></div>
          <div class="br-bar-value">{esc(disp)}</div>
        </div>""")
    return f"""
    <div class="br-card">
      <h2>{esc(title)}</h2>
      <p class="br-card-sub">{esc(subtitle)}</p>
      {"".join(rows)}
      <p class="br-footnote">{esc(footnote)}</p>
    </div>"""


def busco_section(busco):
    if not busco:
        return ""
    segs = [
        ("Complete (single)", busco["S"], "var(--busco-1)"),
        ("Complete (duplicated)", busco["D"], "var(--busco-2)"),
        ("Fragmented", busco["F"], "var(--busco-3)"),
        ("Missing", busco["M"], "var(--busco-4)"),
    ]
    bar = "".join(
        f'<div class="br-stack-seg" style="width:{v}%; background:{c}" title="{esc(l)} — {v}%"></div>'
        for l, v, c in segs if v > 0
    )
    legend = "".join(
        f'<div class="br-legend-item"><span class="br-legend-dot" style="background:{c}"></span>{esc(l)} — {v}%</div>'
        for l, v, c in segs
    )
    n_txt = f", n={fmt(busco['N'])}" if busco.get("N") else ""
    return f"""
    <div class="br-card">
      <h2>Genome BUSCO completeness</h2>
      <p class="br-card-sub">{esc(busco.get('lineage', ''))}{n_txt}</p>
      <div class="br-stack">{bar}</div>
      <div class="br-legend">{legend}</div>
    </div>"""


def kpi_tile(label, value):
    return f"""<div class="br-tile"><div class="br-tile-label">{esc(label)}</div><div class="br-tile-value">{esc(value)}</div></div>"""


def render_pdf(html_path, pdf_path):
    """Print the just-written HTML report to PDF via headless Chromium, so
    the PDF is a faithful rendering of the same file rather than a second,
    separately-maintained layout. Best-effort: warns and skips rather than
    failing the whole report if no Chromium-family browser is on PATH."""
    browser = next(
        (b for b in ("chromium", "chromium-browser", "google-chrome", "google-chrome-stable")
         if shutil.which(b)),
        None,
    )
    if not browser:
        print("Warning: no chromium/google-chrome binary found on PATH — skipping PDF render "
              "(the HTML and .tbl were still written).")
        return False
    abs_html = os.path.abspath(html_path)
    # Chromium needs a writable user-data-dir AND a writable $HOME: running
    # as an arbitrary non-root UID (as Nextflow always does) has no real
    # home directory, so without both it fails outright ("Failed to create
    # headless user data directory container", then a crashpad handler
    # crash even with --user-data-dir alone — crashpad's own database path
    # resolution depends on $HOME regardless). The task's own CWD is always
    # writable (it's the Nextflow work dir), so point both there.
    workdir = os.path.abspath(os.path.dirname(pdf_path) or ".")
    profile_dir = os.path.join(workdir, ".chromium-profile")
    env = {**os.environ, "HOME": workdir}
    result = subprocess.run(
        [browser, "--headless", "--disable-gpu", "--no-sandbox",
         f"--user-data-dir={profile_dir}",
         "--no-pdf-header-footer", f"--print-to-pdf={pdf_path}", f"file://{abs_html}"],
        capture_output=True, text=True, env=env,
    )
    if result.returncode != 0 or not os.path.exists(pdf_path):
        print(f"Warning: PDF render failed (exit {result.returncode}): {result.stderr.strip()[-500:]}")
        return False
    return True


# ----------------------------------------------------------------- main ---
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome-fasta")
    p.add_argument("--final-gff")
    p.add_argument("--proteins-faa")
    p.add_argument("--busco-summary")
    p.add_argument("--annotation-stats")
    p.add_argument("--dbcan-overview")
    p.add_argument("--antismash-json")
    p.add_argument("--phi-base-tsv")
    p.add_argument("--eggnog-annotations")
    p.add_argument("--kegg-module-completeness")
    p.add_argument("--cog-def-tab")
    p.add_argument("--species", required=True)
    p.add_argument("--strain", default="")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    label = f"{args.species}_{args.strain}" if args.strain else args.species
    display_name = f"{args.species.replace('_', ' ')} {args.strain}".strip()

    gstats = genome_stats(args.genome_fasta) if not is_absent(args.genome_fasta) else None
    busco = busco_stats(args.busco_summary) if not is_absent(args.busco_summary) else None
    if not is_absent(args.final_gff):
        feat_counts, gene_lengths, ncrna_lengths = feature_counts(args.final_gff)
    else:
        feat_counts, gene_lengths, ncrna_lengths = None, [], []
    n_proteins = protein_count(args.proteins_faa) if not is_absent(args.proteins_faa) else None
    func_stats = functional_stats(args.annotation_stats) if not is_absent(args.annotation_stats) else None
    cazy_counts, cazy_n = (cazyme_classes(args.dbcan_overview) if not is_absent(args.dbcan_overview) else (None, 0))
    bgc_counts, n_regions = (antismash_bgc_types(args.antismash_json) if not is_absent(args.antismash_json) else (None, 0))
    phi_top, phi_total = (phi_base_phenotypes(args.phi_base_tsv) if not is_absent(args.phi_base_tsv) else (None, 0))
    if not is_absent(args.eggnog_annotations):
        cog_counts, cog_unmapped, cog_none = cog_categories(args.eggnog_annotations, args.cog_def_tab)
    else:
        cog_counts, cog_unmapped, cog_none = None, 0, 0
    kegg_all, kegg_top = (kegg_modules(args.kegg_module_completeness) if not is_absent(args.kegg_module_completeness) else (None, None))

    # ============================================================ .tbl ===
    tbl_path = f"{args.outdir}/{label}_annotation_summary.tbl"
    with open(tbl_path, "w") as f:
        def w(*cols):
            f.write("\t".join(str(c) for c in cols) + "\n")

        w("category", "metric", "value")
        if gstats:
            w("genome", "contigs", gstats["contigs"])
            w("genome", "total_bp", gstats["total_bp"])
            w("genome", "N50", gstats["n50"])
            w("genome", "largest_contig_bp", gstats["largest"])
        if busco:
            w("genome", "busco_lineage", busco.get("lineage", ""))
            w("genome", "busco_complete_single_pct", busco["S"])
            w("genome", "busco_complete_duplicated_pct", busco["D"])
            w("genome", "busco_fragmented_pct", busco["F"])
            w("genome", "busco_missing_pct", busco["M"])
            w("genome", "busco_total_groups", busco.get("N", ""))
        if feat_counts:
            for ftype, n in sorted(feat_counts.items(), key=lambda x: -x[1]):
                w("structural_feature", ftype, n)
            if gene_lengths:
                w("structural_feature", "max_gene_length_bp", max(gene_lengths))
            if ncrna_lengths:
                w("structural_feature", "max_ncRNA_length_bp", max(ncrna_lengths))
        if n_proteins is not None:
            w("structural_feature", "protein", n_proteins)
        if func_stats:
            for lbl, count, pct in func_stats:
                w("functional_annotation", lbl, f"{count} ({pct}%)")
        if cazy_counts:
            for cls, n in sorted(cazy_counts.items(), key=lambda x: -x[1]):
                w("cazyme_class", cls, n)
        if bgc_counts:
            for prod, n in sorted(bgc_counts.items(), key=lambda x: -x[1]):
                w("antismash_bgc_type", prod, n)
            w("antismash_bgc_type", "TOTAL_REGIONS", n_regions)
        if phi_top:
            for pheno, n in phi_top:
                w("phi_base_phenotype", pheno, n)
        if cog_counts:
            for c, n in sorted(cog_counts.items(), key=lambda x: -x[1]):
                w("cog_category", f"{c} ({COG_NAMES[c]})", n)
        if kegg_all:
            for m in kegg_all:
                w("kegg_module_completeness", f"{m['module_id']} {m['module_name']}",
                  f"{m['found']}/{m['total']} ({m['pct']}%)")
    print(f"Wrote {tbl_path}")

    # ============================================================ html ===
    kpis = []
    if gstats:
        kpis.append(kpi_tile("Assembly size", compact_bp(gstats["total_bp"])))
        kpis.append(kpi_tile("N50", compact_bp(gstats["n50"])))
    if busco:
        kpis.append(kpi_tile("BUSCO complete", f"{busco['S'] + busco['D']:.1f}%"))
    if feat_counts:
        kpis.append(kpi_tile("Protein-coding genes", f"{feat_counts.get('gene', 0):,}"))
    if func_stats:
        total_genes = func_stats[0][1]
        with_product = next((c for l, c, pc in func_stats if l == "With product name"), None)
        if with_product is not None and total_genes:
            kpis.append(kpi_tile("With product name", f"{100*with_product/total_genes:.0f}%"))
    if bgc_counts is not None:
        kpis.append(kpi_tile("BGC clusters", f"{n_regions}"))

    struct_items = []
    if feat_counts:
        order = ["gene", "mRNA", "exon", "CDS", "five_prime_UTR", "three_prime_UTR", "ncRNA", "tRNA", "rRNA"]
        struct_items = [(o, feat_counts.get(o, 0), f"{feat_counts.get(o, 0):,}") for o in order if feat_counts.get(o, 0)]

    func_items = []
    if func_stats:
        func_items = [(l, c, f"{pc}% ({c:,})") for l, c, pc in func_stats]

    cazy_items = []
    if cazy_counts:
        cazy_items = [(k, v, f"{v:,}") for k, v in sorted(cazy_counts.items(), key=lambda x: -x[1])]

    bgc_items = []
    if bgc_counts:
        bgc_items = [(k, v, f"{v}") for k, v in sorted(bgc_counts.items(), key=lambda x: -x[1])]

    phi_items = []
    if phi_top:
        phi_items = [(k.replace("__", " + ").replace("_", " "), v, f"{v:,}") for k, v in phi_top]

    cog_items = []
    if cog_counts:
        cog_items = [(f"{c} — {COG_NAMES[c]}", n, f"{n:,}") for c, n in sorted(cog_counts.items(), key=lambda x: -x[1])]

    sections = [
        busco_section(busco),
        bar_section("Structural annotation feature counts", "Feature type breakdown of the final GFF", struct_items),
        bar_section("Functional annotation coverage", f"% of total genes ({fmt(func_stats[0][1]) if func_stats else 0})", func_items),
        bar_section("CAZyme classes", f"dbCAN, {cazy_n:,} annotated genes", cazy_items),
        bar_section("Secondary metabolite BGC types", f"AntiSMASH, {n_regions} regions total", bgc_items),
        bar_section("Top PHI-base pathogenicity phenotypes", f"{phi_total:,} total hits", phi_items),
        bar_section("COG functional categories",
                    f"eggNOG, {sum(cog_counts.values()):,} categorized genes" if cog_counts else "",
                    cog_items,
                    footnote=(f"{cog_unmapped} genes had a raw COG accession not in the reference table and "
                              f"{cog_none} had no COG hit; both excluded above.") if cog_counts else None),
        kegg_module_section(
            "Most complete KEGG pathway modules",
            f"Top {len(kegg_top)} of {len(kegg_all)} annotated modules, min. 4 steps" if kegg_top else "",
            kegg_top,
            "Each module step (one ORTHOLOGY line, e.g. “hexokinase/glucokinase”) counts as satisfied "
            "if the genome has any one of that step's alternative KOs or EC numbers — KEGG catalogs 2–6 "
            "isozyme KOs per step for well-studied pathways, so scoring by raw KO count would penalize exactly "
            "the universal pathways (e.g. glycolysis) most. Still ignores inter-step AND/OR grouping from the "
            "module's full DEFINITION logic, so this remains an approximation, not a formal score.",
        ) if kegg_top else "",
    ]
    sections_html = "".join(s for s in sections if s)

    html_out = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(display_name)} — bagRNA annotation report</title>
<style>{CSS}</style>
</head><body>
<div class="br-root"><div class="br-wrap">
  <div class="br-header">
    <h1>bagRNA Genome Annotation Report</h1>
    <p class="br-species">{esc(display_name)}</p>
    <p class="br-sub">Structural + functional annotation summary</p>
  </div>
  <div class="br-kpis">{"".join(kpis)}</div>
  {sections_html}
  <p class="br-footer">Generated by bagRNA — Nextflow DSL2 genome annotation pipeline</p>
</div></div>
</body></html>"""

    html_path = f"{args.outdir}/{label}_annotation_report.html"
    with open(html_path, "w") as f:
        f.write(html_out)
    print(f"Wrote {html_path}")

    pdf_path = f"{args.outdir}/{label}_annotation_report.pdf"
    if render_pdf(html_path, pdf_path):
        print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
