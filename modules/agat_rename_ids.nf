process AGAT_RENAME_IDS {
    tag "agat_rename_ids"
    label 'process_medium'
    container 'quay.io/biocontainers/agat:1.6.1--pl5321hdfd78af_1'

    input:
    path gff
    val  locus_tag

    output:
    path "agat_renamed.gff3", emit: agat_renamed_gff

    script:
    """
    agat_sp_manage_IDs.pl \\
        --gff ${gff} \\
        --prefix ${locus_tag} \\
        --tair \\
        -o agat_renamed.gff3
    """

    stub:
    """
    touch agat_renamed.gff3
    """
}
