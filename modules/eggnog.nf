process EGGNOG {
    tag "eggnog"
    container 'gabrielerigano/eggnog-mapper:3.0.0-beta6'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa
    val  species
    val  strain
    val  eggnog_data_dir

    output:
    path "eggnog_output.emapper.annotations", emit: eggnog_annotations

    script:
    """
    emapper.py \\
        -i ${proteins_faa} \\
        --itype proteins \\
        -m diamond \\
        --cpu ${task.cpus} \\
        --dmnd_block_size 8 \\
        --data_dir ${eggnog_data_dir} \\
        --output_dir . \\
        -o eggnog_output \\
        --override
    """

    stub:
    """
    touch eggnog_output.emapper.annotations
    """
}
