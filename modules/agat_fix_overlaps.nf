process AGAT_FIX_OVERLAPS {
    tag "${meta.id}"
    label 'process_medium'
    containerOptions '-v $PWD:/data'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    publishDir "${params.outdir}/barrnap", mode: 'copy'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("reformat_barrnap.gff"), emit: barrnap_gff

    script:
    """
    agat_sp_fix_overlaping_genes.pl \\
        --gff /data/${gff} \\
        -o /data/reformat_barrnap.gff
    """

    stub:
    """
    touch reformat_barrnap.gff
    """
}
