process WRAP_GENOME_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("wrapped.fasta"), emit: fasta

    script:
    """
    # table2asn parses '|' in a FASTA header as its own delimited Seq-id
    # syntax (gi|123|gb|ABC, etc.); assembler-generated IDs like Redundans'
    # 'scaffold1|size528146' don't match that shape and get rejected as
    # "not a valid local ID". Sanitize on header lines only, before wrapping.
    seqkit seq -w 60 ${fasta} | sed '/^>/s/|/_/g' > wrapped.fasta
    """

    stub:
    """
    touch wrapped.fasta
    """
}
