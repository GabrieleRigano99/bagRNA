process ANTISMASH {
    tag "antismash"
    label 'process_high'
    container 'antismash/standalone:8.0.4'
    containerOptions '--entrypoint ""'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path gbf_file
    val  species
    val  strain

    output:
    path "antismash_output/", emit: antismash_dir

    script:
    """
    ln -sf ${gbf_file} input.gbk
    antismash \\
        input.gbk \\
        --rre \\
        --genefinding-tool none \\
        -t fungi \\
        --cpus ${task.cpus} \\
        --output-dir antismash_output \\
        --output-basename antismash_${species}_${strain} \\
        --verbose
    """

    stub:
    """
    mkdir -p antismash_output
    """
}
