process ANNEVO_GFF2GTF {
    tag "${meta.id}"
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    publishDir "${params.outdir}/annevo", mode: 'copy'

    input:
    tuple val(meta), path(annevo_gff)

    output:
    tuple val(meta), path("annevo.gtf"), emit: annevo_gtf

    script:
    """
    agat_convert_sp_gff2gtf.pl \\
        --gff ${annevo_gff} \\
        -o annevo.gtf
    """

    stub:
    """
    touch annevo.gtf
    """
}
