process REFORMAT_LOCUS_TAG_IDS {
    tag "reformat_locus_tag_ids"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path agat_renamed_gff
    val  locus_tag
    val  output_name

    output:
    path "${output_name}", emit: renamed_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/reformat_locus_tag_ids.py \\
        ${agat_renamed_gff} \\
        ${locus_tag} \\
        ${output_name}
    """

    stub:
    """
    touch ${output_name}
    """
}
