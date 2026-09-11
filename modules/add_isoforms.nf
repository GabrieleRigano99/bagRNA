process ADD_ISOFORMS {
    tag "add_isoforms"
    label 'process_low'
    container 'python:3.12-slim'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path filtered_pick_gff   // from FILTER_UNSUPPORTED
    path td2_genome_gff3     // TD2 genome-space GFF3 from TD2_GENOME_GFF
    path portcullis_beds     // collected portcullis pass junction BED12 files
    path quant_sf            // Salmon quant.sf for TPM filtering (or NO_FILE to disable)

    output:
    path "isoforms_added.gff3", emit: augmented_gff

    script:
    def min_tpm  = params.add_isoforms_min_tpm != null ? params.add_isoforms_min_tpm : 1.0
    def tpm_flag = quant_sf.name != 'NO_FILE' ? "--quant-sf ${quant_sf} --min-tpm ${min_tpm}" : ""
    """
    python3 ${workflow.projectDir}/bin/add_isoforms.py \\
        ${filtered_pick_gff} \\
        ${td2_genome_gff3} \\
        ${portcullis_beds} \\
        isoforms_added.gff3 \\
        ${tpm_flag}
    """

    stub:
    """
    cp ${filtered_pick_gff} isoforms_added.gff3
    """
}
