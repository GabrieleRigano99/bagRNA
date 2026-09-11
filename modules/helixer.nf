process HELIXER {
    tag "${meta.id}"
    label 'process_high'
    container 'gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1'
    containerOptions "${params.use_gpu ? '--gpus all' : ''}"

    publishDir "${params.outdir}/helixer", mode: 'copy'

    input:
    tuple val(meta), path(fasta)
    val   lineage
    val   species

    output:
    tuple val(meta), path("helixer.gff"), emit: helixer_gff

    script:
    """
    fetch_helixer_models.py -l ${lineage}

    Helixer.py \\
        --fasta-path ${fasta} \\
        --lineage ${lineage} \\
        --gff-output-path helixer.gff \\
        --species ${species}
    """

    stub:
    """
    touch helixer.gff
    """
}
