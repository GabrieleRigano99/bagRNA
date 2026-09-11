process FILTER_UNSUPPORTED {
    tag "filter_unsupported"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path filtered_pick_gff   // from FILTER_ISOFORMS
    path portcullis_beds     // collected portcullis pass junction BED12 files

    output:
    path "unsupported_filtered.gff3", emit: filtered_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/filter_unsupported.py \\
        ${filtered_pick_gff} \\
        ${portcullis_beds} \\
        unsupported_filtered.gff3
    """

    stub:
    """
    cp ${filtered_pick_gff} unsupported_filtered.gff3
    """
}
