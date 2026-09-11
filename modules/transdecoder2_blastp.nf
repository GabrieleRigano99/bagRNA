process TRANSDECODER2_BLASTP {
    tag "td2_blastp"
    label 'process_high'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    input:
    path pep           // longest_orfs.pep from TRANSDECODER2_LONGORFS
    path proteins_dmnd // pre-built protein evidence DIAMOND db (from DIAMOND_MAKEDB)

    output:
    path "td2_blastp_hits.outfmt6", emit: blastp_hits

    script:
    """
    diamond blastp \\
        --query ${pep} \\
        --db ${proteins_dmnd.toString().replaceAll(/\.dmnd$/, '')} \\
        --outfmt 6 \\
        --evalue 1e-5 \\
        --max-target-seqs 1 \\
        --threads ${task.cpus} \\
        --more-sensitive \\
        > td2_blastp_hits.outfmt6
    """

    stub:
    """
    touch td2_blastp_hits.outfmt6
    """
}
