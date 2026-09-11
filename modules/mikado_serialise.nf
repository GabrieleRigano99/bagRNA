process MIKADO_SERIALISE {
    tag "mikado_serialise"
    label 'process_high'
    container 'baderlab/mikado:ubuntu22_mikado2.3.2'

    input:
    path(mikado_dir, stageAs: 'mikado_prepare_dir')
    path(configuration, stageAs: 'configuration.yaml')
    val  scoring_yaml
    path(scoring_file, stageAs: 'scoring_custom.yaml')
    path fasta
    path junctions_bed
    path transdecoder_bed
    path diamond_tsv
    path external_scores
    path prot_evidence

    output:
    path "mikado_op", emit: mikado_dir_serialised

    script:
    """
    cp -rL mikado_prepare_dir mikado_op

    if [ -s "scoring_custom.yaml" ]; then
        cp scoring_custom.yaml ${scoring_yaml}
    else
        cp /usr/local/lib/python3.10/dist-packages/Mikado/configuration/scoring_files/HISTORIC/${scoring_yaml} ./
    fi

    mikado serialise \\
        --json-conf configuration.yaml \\
        --orfs ${transdecoder_bed} \\
        --transcripts mikado_op/mikado_prepared.fasta \\
        --tsv ${diamond_tsv} \\
        --blast-targets ${prot_evidence} \\
        --external-scores ${external_scores} \\
        --genome ${fasta} \\
        --procs ${task.cpus}
    """

    stub:
    """
    cp -rL mikado_prepare_dir mikado_op
    """
}
