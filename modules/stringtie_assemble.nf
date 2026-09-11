process STRINGTIE_ASSEMBLE {
    tag "stringtie_assembly"
    label 'process_high'
    container 'quay.io/biocontainers/stringtie:3.0.3--h29c0135_0'

    publishDir "${params.outdir}/transcript_assembly", mode: 'copy'

    input:
    path bam
    path gtf
    val  strandedness
    path lr_bam  // long-read BAM for --mix mode; NO_FILE to skip

    output:
    path "stringtie.gtf", emit: stringtie_gtf

    script:
    def strand_arg = strandedness == 'secondstrand' ? '--fr' : strandedness == 'firststrand' ? '--rf' : ''
    def mix_arg    = lr_bam.name != 'NO_FILE' ? "--mix ${lr_bam}" : ''
    """
    stringtie \\
        ${bam} \\
        -G ${gtf} \\
        -o stringtie.gtf \\
        -p ${task.cpus} \\
        ${strand_arg} \\
        ${mix_arg} \\
        -j 3 \\
        -v
    """

    stub:
    """
    touch stringtie.gtf
    """
}
