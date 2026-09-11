process STAR_INDEX {
    tag "genome_index"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

    input:
    path fasta
    path gtf

    output:
    path "star_index/", emit: index

    script:
    """
    mkdir -p star_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --genomeSAindexNbases ${params.genomeSAindexNbases} \\
        --runThreadN ${task.cpus}
    """

    stub:
    """
    mkdir -p star_index
    """
}
