process FUNANNOTATE_RENAME {
    tag "funannotate_rename_coding"
    label 'process_medium'
    container 'nextgenusfs/funannotate:v1.8.17'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path coding_gff
    path fasta

    output:
    path "renamed_coding_models_mikado.gff", emit: renamed_gff

    script:
    """
    funannotate gff-rename \\
        -g ${coding_gff} \\
        -o renamed_coding_models_mikado.gff \\
        -f ${fasta}
    """

    stub:
    """
    touch renamed_coding_models_mikado.gff
    """
}
