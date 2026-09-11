process GMAP_ALIGN_TRINITY {
    tag "gmap_trinity"
    label 'process_medium'
    container 'quay.io/biocontainers/gmap:2025.07.31--pl5321hb1d24b7_0'

    publishDir "${params.outdir}/transcript_assembly", mode: 'copy'

    input:
    path transcripts
    path gmap_index

    output:
    path "trinity.gff", emit: trinity_gff

    script:
    """
    gmap.sse42 \\
        -D ${gmap_index} \\
        -d genome \\
        -f gff3_gene \\
        -n 0 \\
        -t ${task.cpus} \\
        --max-intronlength-middle=${params.max_intron_length} \\
        --max-intronlength-ends=${params.max_intron_length} \\
        ${transcripts} \\
        > trinity.gff
    """

    stub:
    """
    touch trinity.gff
    """
}
