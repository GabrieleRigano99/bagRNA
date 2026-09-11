process CPC2 {
    tag "cpc2"
    label 'process_medium'
    container 'ftricomi/cpc2:latest'
    containerOptions '--entrypoint ""'

    publishDir "${params.outdir}/cpc2", mode: 'copy'

    input:
    path fasta

    output:
    path "CPC2_result.txt", emit: cpc2_result

    script:
    """
    python3 \$(which CPC2.py) \\
        -i ${fasta} \\
        -o CPC2_result
    """

    stub:
    """
    touch CPC2_result.txt
    """
}
