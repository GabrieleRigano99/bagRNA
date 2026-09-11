process SAMTOOLS_SPLIT {
    tag "samtools_split"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam)

    output:
    path "split_*.bam", emit: split_bams

    script:
    """
    samtools split \\
        -@ ${task.cpus} \\
        -f 'split_%!.bam' \\
        ${bam}
    """

    stub:
    """
    touch split_sample.bam
    """
}
