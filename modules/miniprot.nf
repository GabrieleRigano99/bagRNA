process MINIPROT {
    tag "miniprot"
    label 'process_high'
    container 'quay.io/biocontainers/compleasm:0.2.8--pyh106432d_0'

    publishDir "${params.outdir}/miniprot", mode: 'copy'

    input:
    path fasta
    path prot_evidence

    output:
    path "miniprot.gtf", emit: miniprot_gtf

    script:
    """
    miniprot ${fasta} ${prot_evidence} \\
        -t ${task.cpus} \\
        --gtf \\
        -G ${params.max_intron_length} \\
        -p 0.70 \\
        -P PEV \\
        > miniprot.gtf
    """

    stub:
    """
    touch miniprot.gtf
    """
}
