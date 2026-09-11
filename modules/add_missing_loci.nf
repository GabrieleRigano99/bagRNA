process ADD_MISSING_LOCI {
    tag "add_missing_loci"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path current_gff       // from GFF_AGAT_FILTER
    path td2_genome_gff3   // from TD2_GENOME_GFF
    path busco_gff         // compleasm genome-mode miniprot hits (COMPLEASM_GENOME)

    output:
    path "struct_augmented.gff3", emit: augmented_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/add_missing_loci.py \\
        ${current_gff} \\
        ${td2_genome_gff3} \\
        ${busco_gff} \\
        struct_augmented.gff3
    """

    stub:
    """
    cp ${current_gff} struct_augmented.gff3
    """
}
