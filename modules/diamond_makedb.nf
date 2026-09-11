process DIAMOND_MAKEDB {
    tag "diamond_makedb"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.1.10--h43eeafb_0'

    publishDir "${params.outdir}/diamond", mode: 'copy'

    input:
    path proteins_fasta

    output:
    path "proteins.dmnd", emit: diamond_db

    script:
    """
    diamond makedb \\
        --in ${proteins_fasta} \\
        -d proteins \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch proteins.dmnd
    """
}
