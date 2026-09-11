process SUBSAMPLE_BAM {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("subsampled.bam"), path("subsampled.bam.bai"), emit: bam_bai

    script:
    """
    samtools view -s 42.05 -b -@ ${task.cpus} ${bam} > subsampled.bam
    samtools index subsampled.bam
    """

    stub:
    """
    touch subsampled.bam subsampled.bam.bai
    """
}
