process DIAMOND_BLASTX {
    tag "diamond_blastx"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.outdir}/mikado", mode: 'copy'

    input:
    path query_fasta
    path diamond_db

    output:
    path "mikado_prepared.diamond.tsv", emit: diamond_tsv

    script:
    """
    diamond blastx \\
        -q ${query_fasta} \\
        -d ${diamond_db} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop \\
        --max-target-seqs 5 \\
        --threads ${task.cpus} \\
        -o mikado_prepared.diamond.tsv
    """

    stub:
    """
    touch mikado_prepared.diamond.tsv
    """
}
