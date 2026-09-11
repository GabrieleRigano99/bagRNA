process TRINITY {
    tag "trinity_assembly"
    label 'process_long'
    container 'trinityrnaseq/trinityrnaseq:2.15.2'

    publishDir "${params.outdir}/trinity", mode: 'copy'

    input:
    path   bam
    val    jaccard_clip
    val    ram_limit
    val    orientation

    output:
    path "trinity_GG.Trinity-GG.fasta", emit: trinity_fasta

    script:
    def jaccard_arg  = jaccard_clip ? '--jaccard_clip' : ''
    def lib_type_arg = orientation  ? "--SS_lib_type ${orientation}" : ''
    """
    Trinity \\
        --genome_guided_bam ${bam} \\
        --genome_guided_max_intron ${params.max_intron_length} \\
        --max_memory ${ram_limit} \\
        --min_contig_length 200 \\
        --normalize_reads \\
        --normalize_max_read_cov ${params.trinity_max_cov} \\
        --min_kmer_cov 2 \\
        ${lib_type_arg} \\
        ${jaccard_arg} \\
        --CPU ${task.cpus} \\
        --full_cleanup \\
        --grid_node_CPU ${task.cpus} \\
        --output trinity_GG
    """

    stub:
    """
    touch trinity_GG.Trinity-GG.fasta
    """
}
