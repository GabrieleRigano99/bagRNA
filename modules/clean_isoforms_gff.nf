process CLEAN_ISOFORMS_GFF {
    tag "clean_isoforms_gff"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path augmented_gff   // from ADD_ISOFORMS

    output:
    path "isoforms_cleaned.gff3", emit: cleaned_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/clean_isoforms_gff.py \\
        ${augmented_gff} \\
        isoforms_cleaned.gff3
    """

    stub:
    """
    cp ${augmented_gff} isoforms_cleaned.gff3
    """
}
