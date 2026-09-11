process EFFECTORP3 {
    tag "effectorp3"
    label 'process_medium'
    container 'gabrielerigano/effectorp3:latest'
    containerOptions '--entrypoint ""'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "effectorp3_output.txt", emit: effectorp3_results

    script:
    """
    python /opt/EffectorP-3.0/EffectorP.py \\
        -f \\
        -i ${proteins_faa} \\
        -o effectorp3_output.txt
    """

    stub:
    """
    touch effectorp3_output.txt
    """
}
