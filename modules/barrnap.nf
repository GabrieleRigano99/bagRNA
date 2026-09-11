process BARRNAP {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/barrnap:0.9--1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("barrnap_raw.gff"), emit: barrnap_gff

    script:
    """
    barrnap \\
        --kingdom euk \\
        --reject 0.50 \\
        --threads ${task.cpus} \\
        ${fasta} \\
    | awk 'BEGIN{FS=OFS="\\t"} /^#/ || \$9 !~ "partial"' \\
    > barrnap_raw.gff
    """

    stub:
    """
    touch barrnap_raw.gff
    """
}
