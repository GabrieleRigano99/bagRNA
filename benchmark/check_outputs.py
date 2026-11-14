#!/usr/bin/env python3
"""Check finished benchmark runs against the pass rule. Standard library only.

Usage: benchmark/check_outputs.py RUNS_DIR [species ...] [--busco BUSCO_DIR]

RUNS_DIR holds <species>/results (the RUNS path of run_species.sh). BUSCO_DIR, if given, holds
<species>/summary.txt written by `compleasm protein -o <species>` on structural_annotation/final_proteins.faa.

A species passes when:
  1. every task in pipeline_info/execution_trace.txt is COMPLETED
  2. the table2asn .val file has 0 "Error:" lines (Warnings are allowed)
  3. gene count is within 10 % of the baseline
  4. % of genes with a KEGG KO and with an InterPro hit are within 2 points of the baseline
  5. (with --busco) complete protein BUSCO is not more than 1 point BELOW the baseline
     (a higher score is reported but does not fail)
Exit code 0 when every checked species passes.
"""
import csv, glob, os, re, sys

# genes (annotated.gff3), %KO, %InterPro, protein BUSCO complete: first bagRNA run (Sep 2026) of each species
BASE = {
    "Aspergillus_aculeatinus":   (17434, 39.3, 57.2, 98.59),
    "Candida_orthopsilosis":     (6358,  67.5, 79.3, 99.04),
    "Cryptococcus_depauperatus": (8192,  52.4, 63.0, 94.42),
    "Malassezia_globosa":        (4904,  70.2, 80.2, 82.40),
    "Sporothrix_schenckii":      (17123, 31.7, 44.0, 96.08),
    "Syncephalastrum_racemosum": (13900, 51.0, 65.3, 97.19),
}


def pct_with(tsv, col):
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    return 100 * sum(r[col] not in ("", "-", "NA") for r in rows) / len(rows)


def busco_complete(summary):
    m = {k: float(v) for k, v in re.findall(r"^([SD]):([\d.]+)%", open(summary).read(), re.M)}
    return m["S"] + m["D"]


def check(sp, runs, busco_dir):
    res = os.path.join(runs, sp, "results")
    out, ok = [], True

    def rec(name, good, text):
        nonlocal ok
        ok &= good
        out.append(f"  {'ok  ' if good else 'FAIL'} {name}: {text}")

    trace = os.path.join(res, "pipeline_info", "execution_trace.txt")
    status = [l.split("\t")[4] for l in open(trace).read().splitlines()[1:]]
    bad = sum(s != "COMPLETED" for s in status)
    rec("tasks", bad == 0, f"{len(status)} tasks, {bad} not COMPLETED")

    val = glob.glob(os.path.join(res, "ncbi_submission", "*.val"))
    text = open(val[0]).read() if val else ""
    n_err, n_warn = text.count("Error:"), text.count("Warning:")
    rec("table2asn", bool(val) and n_err == 0, f"{n_err} Error, {n_warn} Warning")

    base_genes, base_ko, base_ipr, base_busco = BASE[sp]
    fa = os.path.join(res, "functional_annotation")
    genes = sum(1 for l in open(os.path.join(fa, "annotated.gff3")) if l.split("\t")[2:3] == ["gene"])
    rec("genes", abs(genes - base_genes) <= 0.10 * base_genes,
        f"{genes} (baseline {base_genes}, {100 * (genes - base_genes) / base_genes:+.1f} %)")
    tsv = os.path.join(fa, "functional_annotation.tsv")
    for name, col, base in (("%KO", "KEGG_KO", base_ko), ("%InterPro", "InterPro_accessions", base_ipr)):
        v = pct_with(tsv, col)
        rec(name, abs(v - base) <= 2, f"{v:.1f} (baseline {base})")

    if busco_dir and not os.path.exists(os.path.join(busco_dir, sp, "summary.txt")):
        out.append("  skip BUSCO protein: no summary.txt in " + busco_dir)
    elif busco_dir:
        b = busco_complete(os.path.join(busco_dir, sp, "summary.txt"))
        rec("BUSCO protein", b >= base_busco - 1, f"{b:.2f} (baseline {base_busco}, {b - base_busco:+.2f})")
    print(f"{sp}: {'PASS' if ok else 'FAIL'}")
    print("\n".join(out))
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    busco = args.pop(args.index("--busco") + 1) if "--busco" in args else None
    if busco:
        args.remove("--busco")
    if not args:
        sys.exit(__doc__)
    runs, species = args[0], args[1:] or list(BASE)
    sys.exit(0 if all([check(sp, runs, busco) for sp in species]) else 1)
