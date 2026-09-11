process PHI_BASE {
    tag "phi_base"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.outdir}/functional_annotation/phi_base", mode: 'copy'

    input:
    path proteins
    path phi_base_db   // FASTA or pre-built .dmnd

    output:
    path "phi_base.tsv", emit: phi_base_tsv

    script:
    def is_dmnd = phi_base_db.name.endsWith('.dmnd')
    def db_arg  = is_dmnd ? phi_base_db.toString().replaceAll(/\.dmnd$/, '') : 'phi_base_db'
    def build   = is_dmnd ? '' : "diamond makedb --in ${phi_base_db} -d phi_base_db --threads ${task.cpus}"
    """
    ${build}
    diamond blastp \\
        --query ${proteins} \\
        --db ${db_arg} \\
        --out phi_base.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \\
        --evalue 1e-10 \\
        --max-target-seqs 1 \\
        --more-sensitive \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch phi_base.tsv
    """
}
