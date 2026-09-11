process MEROPS {
    tag "merops"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.outdir}/functional_annotation/merops", mode: 'copy'

    input:
    path proteins
    path merops_db   // FASTA or pre-built .dmnd

    output:
    path "merops.tsv", emit: merops_tsv

    script:
    def is_dmnd = merops_db.name.endsWith('.dmnd')
    def db_arg  = is_dmnd ? merops_db.toString().replaceAll(/\.dmnd$/, '') : 'merops_db'
    def build   = is_dmnd ? '' : "diamond makedb --in ${merops_db} -d merops_db --threads ${task.cpus}"
    """
    ${build}
    diamond blastp \\
        --query ${proteins} \\
        --db ${db_arg} \\
        --out merops.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \\
        --evalue 1e-5 \\
        --max-target-seqs 1 \\
        --more-sensitive \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch merops.tsv
    """
}
