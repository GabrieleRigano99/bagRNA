process SALMON_QUANT {
    tag "salmon_quant"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'

    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    path salmon_index
    path reads_r1_list
    path reads_r2_list

    output:
    path "salmon_quant/quant.sf", emit: quant_sf

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads_r1_list} \\
        -2 ${reads_r2_list} \\
        -p ${task.cpus} \\
        --seqBias \\
        --posBias \\
        --gcBias \\
        --validateMappings \\
        -o salmon_quant
    """

    stub:
    """
    mkdir -p salmon_quant
    touch salmon_quant/quant.sf
    """
}
