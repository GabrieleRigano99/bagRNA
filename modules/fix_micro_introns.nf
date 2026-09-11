process FIX_MICRO_INTRONS {
    tag "fix_micro_introns"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path recovered_gff   // from ADD_BUSCO_ISOFORMS
    path fasta

    output:
    path "struct_microfixed.gff3", emit: fixed_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/fix_micro_introns.py \\
        ${recovered_gff} \\
        ${fasta} \\
        struct_microfixed.gff3
    """

    stub:
    """
    cp ${recovered_gff} struct_microfixed.gff3
    """
}
