process MIKADO_PICK {
    tag "mikado_pick"
    label 'process_high'
    container 'baderlab/mikado:ubuntu22_mikado2.3.2'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path(mikado_dir, stageAs: 'mikado_op_input')
    path(configuration, stageAs: 'configuration.yaml')
    val  scoring_yaml
    path(scoring_file, stageAs: 'scoring_custom.yaml')
    path fasta

    output:
    path "mikado_op/mikado_pick.loci.gff3",    emit: pick_gff
    path "mikado_op/mikado_monoloci.gff3",      emit: monoloci_gff
    path "mikado_op/mikado_subloci.gff3",       emit: subloci_gff

    script:
    """
    cp -rL mikado_op_input mikado_op

    if [ -s "scoring_custom.yaml" ]; then
        cp scoring_custom.yaml ${scoring_yaml}
    else
        cp /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/HISTORIC/${scoring_yaml} ./
    fi

    mikado pick \\
        --json-conf configuration.yaml \\
        --loci-out mikado_pick.loci.gff3 \\
        --subloci-out mikado_subloci.gff3 \\
        --monoloci-out mikado_monoloci.gff3 \\
        --genome ${fasta} \\
        --procs ${task.cpus}
    """

    stub:
    """
    cp -rL mikado_op_input mikado_op
    touch mikado_op/mikado_pick.loci.gff3
    touch mikado_op/mikado_subloci.gff3
    touch mikado_op/mikado_monoloci.gff3
    """
}
