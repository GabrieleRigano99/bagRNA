process ADD_BUSCO_ISOFORMS {
    tag "add_busco_isoforms"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy'

    input:
    path augmented_gff   // from ADD_MISSING_LOCI
    path busco_gff       // compleasm genome-mode miniprot hits (COMPLEASM_GENOME)
    path missing_ids     // from COMPLEASM_MISSING_BUSCO

    output:
    path "struct_recovered.gff3", emit: recovered_gff

    script:
    """
    python3 ${workflow.projectDir}/bin/add_busco_isoforms.py \\
        ${augmented_gff} \\
        ${busco_gff} \\
        ${missing_ids} \\
        struct_recovered.gff3
    """

    stub:
    """
    cp ${augmented_gff} struct_recovered.gff3
    """
}
