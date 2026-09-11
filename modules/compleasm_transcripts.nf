process COMPLEASM_TRANSCRIPTS {
    tag "busco_protein_scores"
    label 'process_high'
    container 'quay.io/biocontainers/compleasm:0.2.8--pyh106432d_0'

    publishDir "${params.outdir}/mikado", mode: 'copy', pattern: "busco_scores.tsv"

    input:
    path td2_pep
    val  lineage
    path busco_db, stageAs: 'mb_downloads'

    output:
    path "busco_scores.tsv", emit: busco_scores

    script:
    def lineage_full = "${lineage}_odb12"
    """
    compleasm protein \\
        -p ${td2_pep} \\
        -l ${lineage_full} \\
        -L mb_downloads \\
        -o compleasm_prot_out \\
        -t ${task.cpus}

    # compleasm protein mode writes full_table.tsv flat under the output dir
    # (not nested under a lineage subdir as genome mode does); locate it robustly.
    full_table=\$(find compleasm_prot_out -name full_table.tsv | head -1)

    python3 ${workflow.projectDir}/bin/busco_scores.py \\
        "\$full_table" \\
        mb_downloads/${lineage_full}/scores_cutoff \\
        > busco_scores.tsv
    """

    stub:
    """
    echo -e "tid\\tbusco_score" > busco_scores.tsv
    """
}
