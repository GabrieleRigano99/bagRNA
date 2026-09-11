process COMPLEASM_MISSING_BUSCO {
    tag "missing_busco_check"
    label 'process_high'
    container 'quay.io/biocontainers/compleasm:0.2.8--pyh106432d_0'

    publishDir "${params.outdir}/structural_annotation", mode: 'copy', pattern: "missing_busco_ids.txt"

    input:
    path interim_proteins
    val  lineage
    path busco_db, stageAs: 'mb_downloads'

    output:
    path "missing_busco_ids.txt", emit: missing_ids

    script:
    def lineage_full = "${lineage}_odb12"
    """
    compleasm protein \\
        -p ${interim_proteins} \\
        -l ${lineage_full} \\
        -L mb_downloads \\
        -o compleasm_prot_out \\
        -t ${task.cpus}

    full_table=\$(find compleasm_prot_out -name full_table.tsv | head -1)
    awk -F'\\t' '\$2=="Missing"{print \$1}' "\$full_table" > missing_busco_ids.txt
    """

    stub:
    """
    touch missing_busco_ids.txt
    """
}
