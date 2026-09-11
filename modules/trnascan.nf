process TRNASCAN {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/trnascan-se:2.0.12--pl5321h7b50bb2_2'

    publishDir "${params.outdir}/trnascan", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("trnascan.gff"), emit: trnascan_gff

    script:
    """
    tRNAscan-SE -E \\
        ${fasta} \\
        --gff trnascan.gff \\
        --thread ${task.cpus}
    """

    stub:
    """
    touch trnascan.gff
    """
}
