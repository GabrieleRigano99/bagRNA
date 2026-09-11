process KEGG_ANNOTATE {
    tag "kegg_annotate"
    label 'process_low'
    container 'python:3.11-slim'
    containerOptions '-v /usr/bin/ps:/usr/bin/ps'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path kofamscan_tsv

    output:
    path "kegg_annotations.tsv", emit: kegg_tsv

    script:
    """
    kegg_annotate.py ${kofamscan_tsv} kegg_annotations.tsv
    """

    stub:
    """
    touch kegg_annotations.tsv
    """
}
