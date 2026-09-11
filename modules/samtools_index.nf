process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: bam_bai

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """

    stub:
    """
    touch ${bam}.bai
    """
}
