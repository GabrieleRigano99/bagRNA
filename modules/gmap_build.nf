process GMAP_BUILD {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/gmap:2025.07.31--pl5321hb1d24b7_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("gmap_index/"), emit: gmap_index

    script:
    """
    gmap_build \\
        -D gmap_index \\
        -d genome \\
        -k 13 \\
        ${fasta}
    """

    stub:
    """
    mkdir -p gmap_index
    """
}
