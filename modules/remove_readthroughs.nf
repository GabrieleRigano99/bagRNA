process REMOVE_READTHROUGHS {
    tag "remove_readthroughs"
    label 'process_low'
    container 'python:3.11-slim'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path mikado_pick_gff
    path miniprot_gtf

    output:
    path "readthrough_filtered.gff3", emit: filtered_gff
    path "readthroughs.txt",          emit: report

    script:
    """
    remove_readthroughs.py \\
        --mikado   ${mikado_pick_gff} \\
        --miniprot ${miniprot_gtf} \\
        --out      readthrough_filtered.gff3 \\
        --report   readthroughs.txt
    """

    stub:
    """
    touch readthrough_filtered.gff3 readthroughs.txt
    """
}
