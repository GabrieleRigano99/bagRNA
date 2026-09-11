process SALMON_INDEX {
    tag "salmon_index"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0'

    input:
    path fasta

    output:
    path "salmon_index/", emit: salmon_index

    script:
    """
    # Salmon ≥1.10 crashes on mixed-case (soft-masked) FASTA; uppercase first
    awk '/^>/{print; next}{print toupper(\$0)}' ${fasta} > fasta_upper.fa
    salmon index \\
        -t fasta_upper.fa \\
        -i salmon_index \\
        -p ${task.cpus}
    """

    stub:
    """
    mkdir -p salmon_index
    """
}
