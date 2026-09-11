process MINIMAP2_ALIGN {
    tag "lr_align"
    label 'process_high'
    container 'rkimf1/minimap2-samtools:2.28--4176384'

    publishDir "${params.outdir}/long_reads", mode: 'copy'

    input:
    path fasta
    path reads         // one or more FASTQ files (collected)
    path junctions_bed // portcullis BED for splice guidance; NO_FILE to skip

    output:
    path "lr_sorted.bam",     emit: bam
    path "lr_sorted.bam.bai", emit: bai

    script:
    def preset   = params.lr_type == 'pacbio_hifi' ? 'splice:hq' : 'splice'
    def junc_arg = junctions_bed.name != 'NO_FILE' ? "--junc-bed ${junctions_bed} --junc-bonus 15" : ''
    """
    minimap2 \\
        -ax ${preset} \\
        --secondary=no \\
        -G ${params.max_intron_length} \\
        ${junc_arg} \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -o lr_sorted.bam

    samtools index lr_sorted.bam
    """

    stub:
    """
    touch lr_sorted.bam lr_sorted.bam.bai
    """
}
