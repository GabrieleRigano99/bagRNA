process GMAP_ALIGN_TRANSCRIPTS {
    tag "gmap_transcript_evidence"
    label 'process_medium'
    container 'quay.io/biocontainers/gmap:2025.07.31--pl5321hb1d24b7_0'

    publishDir "${params.outdir}/transcript_evidence", mode: 'copy'

    input:
    path transcripts
    path gmap_index

    output:
    path "transcript_evidence.gff", emit: transcript_evidence_gff

    script:
    """
    gmap.sse42 \\
        -D ${gmap_index} \\
        -d genome \\
        -f gff3_gene \\
        -n 0 \\
        -t ${task.cpus} \\
        ${transcripts} \\
        > transcript_evidence.gff
    """

    stub:
    """
    touch transcript_evidence.gff
    """
}
