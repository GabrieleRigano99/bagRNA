process DEEPKOALA {
    tag "deepkoala"
    label 'process_medium'
    container 'gabrielerigano/deepkoala:0.1-beta'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "kofamscan_result.tsv", emit: kofamscan_tsv

    script:
    """
    deepkoala \\
        -i ${proteins_faa} \\
        -o deepkoala_raw.csv \\
        --model full \\
        --detail \\
        --device auto \\
        --batch_size 128 \\
        --num_workers ${task.cpus}

    # Reformat into KofamScan's own detail-tsv column layout so every
    # downstream consumer (KEGG_ANNOTATE, KEGG_KO_INFO, ANNOTATE_FUNCTIONAL's
    # merge_functional_annotations.py) works unmodified: they only ever read
    # column 0 (significance marker), 1 (gene_id) and 2 (KO). Column 6
    # (description) is deliberately left EMPTY rather than a placeholder —
    # merge_functional_annotations.py only overwrites the product/EC fields
    # from this file when that column is non-empty; leaving it blank lets
    # the later KEGG REST enrichment step (fetch_kegg_ko_info.py, keyed on
    # the same KO) fill in product name and EC number from KEGG's own NAME
    # field instead, which DeepKOALA has no equivalent of natively.
    python3 - <<'PYEOF'
import csv
with open("deepkoala_raw.csv") as fin, open("kofamscan_result.tsv", "w") as fout:
    for row in csv.DictReader(fin):
        if row["annotate"] != "*":
            continue
        fout.write("\\t".join([
            "*", row["name"], row["predict_label"],
            row["threshold"], row["probability"], "-", "",
        ]) + "\\n")
PYEOF
    """

    stub:
    """
    touch kofamscan_result.tsv
    """
}
