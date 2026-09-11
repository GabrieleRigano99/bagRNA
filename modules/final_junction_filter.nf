process FINAL_JUNCTION_FILTER {
    tag "final_junction_filter"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path cleaned_gff       // from CLEAN_ISOFORMS_GFF
    path portcullis_beds   // collected portcullis pass junction BED12 files

    output:
    path "isoforms_juncfiltered.gff3", emit: filtered_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/final_junction_filter.py \\
        ${cleaned_gff} \\
        ${portcullis_beds} \\
        isoforms_juncfiltered.gff3
    """

    stub:
    """
    cp ${cleaned_gff} isoforms_juncfiltered.gff3
    """
}
